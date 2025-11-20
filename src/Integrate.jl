module Integrate
    using Distributed
    using OrdinaryDiffEq
    using FLoops
    using DistributedArrays
    include("math.jl")
    include("save.jl")
    
    """
        Returns initial conditions for the beginning of Roche Lobe overflow as a vector.
        Vector output is [a, M_WD, M_NS, R_WD, J, θ], where:
        a = total orbital separation of the binary.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        R_WD = radius of the white dwarf.
        J = total angular momentum of the binary.
        θ = angular parameter.
        ```julia
        initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2.)
        ```
        mu_e = mean molecular mass per electron (default 2 for white dwarfs of all but the most exotic chemical compositions (e.g., ONeMG)).
    """
    function initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2., mode::String="Debug")::Vector{Float64}
        #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
        #M_NS=1.35 M_\odot is typical NS mass
        if mode=="Debug"
            @assert M_WD<1.4 #Chandrasekar mass, beyond which we cannot calculate R_WD(see initial function)
        end
        a=Math.RL_contact(M_WD,M_WD,mu_e=mu_e)
        J=Math.circular_J(a,M_WD,M_WD)
        R_WD=Math.R0_WD(M_WD)
        return([a, M_WD, M_NS, R_WD, J, 0.])
    end
    
    #"Below is the heart of the integration process"
    """Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken and the updated theta parameters."""
    function one_step!(theta::Vector{Float64}; dt_period::Float64=1e-2, dt::Float64=-1., mode::String="Debug", parallel::Bool=true, uneven::Bool=false)::Tuple{Float64, Vector{Float64}}
        if mode=="Debug"
            @assert length(theta)==6
        end
        #integrates forward one step with time interval either given explicitly (dt!=-1)
        #or as a fraction of the orbital period determined by dt_period
        if dt==-1
            @inbounds begin
                a = theta[1]
                M_WD = theta[2]
                M_NS = theta[3]
                dt=dt_period*Math.period(a, M_WD + M_NS)
            end
        end

        prob = nothing
        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        redirect_stderr(devnull) do 
            if uneven
                prob = ODEProblem(Math.acceleration_uneven!,theta,(0.,dt),[])
            elseif parallel
                prob = ODEProblem(Math.acceleration_even!,theta,(0.,dt),[])
            else
                prob = ODEProblem(Math.acceleration_serial!,theta,(0.,dt),[])
            end
        end
        theta .= solve(prob, alg=Vern9()).u[end]
        dt,theta
    end
    
    """Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. The function also saves the run to a file named by the date and initial conditions if `write_file`. """
    function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3, dt_period=1e-3,
            parallel::Bool=true, uneven::Bool=false, write_file::Bool=true, filename::Union{String,Int}=-1)::AbstractMatrix{<:Real}
        #returns time series second derivatives of the quadrupole moment
        @assert M_WD > mass_limit
        @assert M_WD < Math.M_ch # Chandrasekhar mass
        @assert M_NS > M_WD
        @assert M_NS > 1.1 # very low lower limit, slightly below mass of the secondary neutron star in PSR J0453+1559
        @assert M_NS < 2.9 # TOV limit (a generous estimate)
        @assert orbit_limit > 0
        @assert step_limit > 1
        @assert mass_limit > 0
        @assert dt_period > 0
        t=0.
        theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
        ddI=Math.ddI_from_theta(theta)
        datumses=zeros(step_limit+1,9)
        datumses[1,:].=[t, ddI[1,1], ddI[1,2], theta... ]
        i=1

        @inbounds begin
            phase = theta[6]
            M_WD = theta[2]
            while phase<2*π*orbit_limit && i<step_limit && M_WD>mass_limit
                #integrate until time or step limit are elapsed or WD completely dissipates
                
                dt,theta=one_step!(theta; dt_period=dt_period, parallel=parallel, uneven=uneven)
                ddI=Math.ddI_from_theta(theta)
                i+=1
                t+=dt
            
                datumses[i,:].=[t, ddI[1,1], ddI[1,2], theta... ]
            end
        end
    
        datumses=datumses[1:i,:]
    
        if write_file #saves datumses to a csv in directory "Output"
            Save.save_as_csv(datumses,M_WD,M_NS,filename=filename)
        end
        
        return(datumses)
    end

    """Performs many WD-NS merger simulations at once, with user specified initial white dwarf and neutron star masses. (This is another place where we take advantage of parallelization). Takes in sim_conditions (NamedTuple), which has:
    - M_WD_lower (lower limit for white dwarf mass)
    - M_WD_upper (upper limit for white dwarf mass)
    - M_WD_step (step size for white dwarf mass)
    - M_NS_lower (lower limit for neutron star mass)
    - M_NS_upper (upper limit for neutron star mass)
    - M_NS_step (step size for neutron star mass)
    The ranges of M_WD and M_NS must be physical. Some of the checks we have are:
    @assert M_WD < Math.M_ch # Chandrasekhar mass
    @assert M_NS > M_WD
    @assert M_NS > 1.1 # very low lower limit, slightly below mass of the secondary neutron star in PSR J0453+1559
    @assert M_NS < 2.9 # TOV limit (a generous estimate)"""
    function many_integrate(sim_conditions::NamedTuple; mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3, dt_period=1e-3,
            parallel::Bool=true, distributed::Bool=true, uneven::Bool=false, write_file::Bool=true, filename::Union{String,Int}=-1)
        @assert sim_conditions.M_WD_lower < sim_conditions.M_WD_upper
        @assert sim_conditions.M_NS_lower < sim_conditions.M_NS_upper
        @assert sim_conditions.M_WD_step <= sim_conditions.M_WD_upper - sim_conditions.M_WD_lower
        @assert sim_conditions.M_NS_step <= sim_conditions.M_NS_upper - sim_conditions.M_NS_lower
        M_WDs = collect(range(sim_conditions.M_WD_lower, stop=sim_conditions.M_WD_upper, step=sim_conditions.M_WD_step))
        M_NSs = collect(range(sim_conditions.M_NS_lower, stop=sim_conditions.M_NS_upper, step=sim_conditions.M_NS_step))
        if parallel
            if distributed
                pairs = [(M_WD, M_NS) for M_WD in M_WDs for M_NS in M_NSs if M_WD < M_NS]
                pairs_distr = distribute(pairs)
                @distributed for pair in pairs_distr
                    M_WD, M_NS = pair
                    integrate(M_WD; M_NS, mu_e, 
                        orbit_limit, step_limit, mass_limit, dt_period,
                        parallel, uneven, write_file, filename)
                    GC.gc() # Free up memory
                end    
                # Need to be careful here because apparently having each process read and write to disk can worsen performance if not done carefully
            else
                @floop ThreadedEx(basesize = 8) for M_WD in M_WDs, M_NS in M_NSs
                    if M_WD < M_NS
                        integrate(M_WD; M_NS, mu_e, 
                            orbit_limit, step_limit, mass_limit, dt_period,
                            parallel, uneven, write_file, filename)
                    end
                    GC.gc() # Free up memory
                end
            end
        else
            for M_WD in M_WDs, M_NS in M_NSs
                if M_WD < M_NS
                    integrate(M_WD; M_NS, mu_e, 
                        orbit_limit, step_limit, mass_limit, dt_period,
                        parallel, uneven, write_file, filename)
                end
                GC.gc() # Free up memory
            end
        end
    end
end