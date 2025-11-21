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


    const M_ch=1.44 #Chandresekhar mass in M⊙
    optional_default=[2.0,1,10.0,false,1.0,"Isotropic Re-emission",0.1,true]
    #`optional` contains optional parameters. In order, these optional parameters are:
    #mu_e=2.0:                     the chemical composition of the white dwarf (for CO/ONeMg, mu_e=2.0/2.35, respectively)
    #equation=1:                   choice of empirical equation for white dwarf mass vs. radius relation (1=Nauenberg/2=Eggleton
    #A=10.0:                       the rate at which mass flows off the WD surface (characteristic number of orbits for mass to leave)
    #poly_shell=true:              if true/false, uses polytropic/uniform density distribution to model mass in process of being stripped from WD
    #a_ring_frac=1.0:              radius of circumbinary ring formed by matter decreted by WD
    #mode="Isotropic Re-emission": mode of mass loss assumed in calculating angular momentum loss
    #eta_acc=0.1:                  accretion efficiency of the NS
    #period_calc=true:             if true/false, calculates angular velocity as a function of radius/angular momentum

    #Of these, I expect A, poly_shell, mode, and period_calc to have considerable effect on the evolution of the system

    function initial_circular_orbit(M_WD::Float64; M_NS::Float64=1.35, optional=optional_default, mode::String="Debug")::Vector{Float64}
    
        #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
        #M_NS=1.35 M_\odot is typical NS mass
        if mode=="Debug"
            @assert M_WD<M_ch #Chandrasekar mass, beyond which we cannot calculate R_WD(see initial function)
        end
        a=Math.RL_contact(M_WD,M_WD;optional=optional)
        J=Math.circular_J(a,M_WD,M_WD)
        R_WD=Math.R0_WD(M_WD)
        return([a, M_WD, M_NS, R_WD, J, 0.])
    end
    
    #"Below is the heart of the integration process"
    """Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken and the updated theta parameters."""
    function one_step!(theta::Vector{Float64}; optional=optional_default, dt_period::Float64=1e-2, dt::Float64=-1., mode::String="Debug", parallel::Bool=true, uneven::Bool=false)::Tuple{Float64, Vector{Float64}}
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
            elseif parallel #define a function with fixed optional parameters
                function t_acceleration!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64)
                    Math.acceleration_even!(du, theta, p, t; optional=optional)
                end
                prob = ODEProblem(t_acceleration!,theta,(0.,dt),[])
            else
                function t_acceleration!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64)
                    Math.acceleration_serial!(du, theta, p, t; optional=optional)
                end
                prob = ODEProblem(t_acceleration!,theta,(0.,dt),[])
            end
        end
        theta .= solve(prob, alg=Vern9()).u[end]
        dt,theta
    end
    
    """Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. The function also saves the run to a file named by the date and initial conditions if `write_file`. """
    function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3, dt_period=1e-3,
            parallel::Bool=true, uneven::Bool=false, write_file::Bool=true, filename::Union{String,Int}=-1,
        optional=optional_default)::AbstractMatrix{<:Real}
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
        theta=initial_circular_orbit(M_WD,M_NS=M_NS,optional=optional)
        ddI=Math.ddI_from_theta(theta,optional=optional)
        datumses=zeros(step_limit+1,9)
        datumses[1,:].=[t, ddI[1,1], ddI[1,2], theta... ]
        i=1

        @inbounds begin
            phase = theta[6]
            M_WD = theta[2]
            while phase<2*π*orbit_limit && i<step_limit && M_WD>mass_limit
                #integrate until time or step limit are elapsed or WD completely dissipates
                
                dt,theta=one_step!(theta; dt_period=dt_period, parallel=parallel, uneven=uneven, optional=optional)
                ddI=Math.ddI_from_theta(theta; optional=optional)
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

    """Performs many WD-NS merger simulations at once. (This is another place where we take advantage of parallelization). Takes in sim_conditions (NamedTuple), which has:
    - M_WDs (array/tuple of white dwarf mass(es))
    - M_NSs (array/tuple of neutron star mass(es))
    and optional_conditions (NamedTuple), which has:
    - mu_es (array/tuple of mu_e(s), the mean molecular mass per electron)
    - WD_mass_radius_equations ((1), (2), or (1,2); specifies whether to use the Nauenberg WD mass-radius relation (1), Eggleton relation (2), or both (1,2))
    - As (array/tuple of A value(s), a constant related to mass loss of the white dwarf)
    - use_polytropic_density ((true), (false), or (true,false); specifies whether to use the polytropic density distribution to model mass in process of being stripped from WD (true), the uniform density distribution (false), or both (true,false))
    - a_ring_fracs (array/tuple of a_ring_frac(s), the radius of circumbinary ring formed by matter decreted by WD. Only relevant if 'modes' includes "Circumbinary Ring")
    - mass_loss_modes (array/tuple of mass loss mode specifications for calculating angular momentum loss; possible values include 'j', 'i', and 'c'. 'j' = Jeans mode, 'i' = isotropic re-emission, 'c' = circumbinary ring)
    - eta_accs (array/tuple of eta_acc(s), the accretion efficiency of the neutron star)
    - use_period_calc ((true), (false), or (true,false); specifies whether to calculate the angular velocity of white dwarf as a function of the period (true), the angular momentum (false), or both (true,false))
    many_integrate then iterates through all possible combinations of these user-specified parameters, with one integrate run dedicated to each combination. Any parameters not specified in optional_conditions will fall back to the default values, which are:
    - mu_es = (2.0)
    - WD_mass_radius_equations = (1)
    - As = (10.0)
    - use_polytropic_density = (true)
    - a_ring_fracs = (1.0)
    - mass_loss_modes = ("Isotropic Re-emission")
    - eta_accs = (0.1)
    - use_period_calc = (true)
    The ranges of M_WD and M_NS must be physical. Some of the checks we have are:
    @assert M_WD < Math.M_ch # Chandrasekhar mass
    @assert M_NS > M_WD
    @assert M_NS > 1.1 # very low lower limit, slightly below mass of the secondary neutron star in PSR J0453+1559
    @assert M_NS < 2.9 # TOV limit (a generous estimate)"""
    function many_integrate(sim_conditions::NamedTuple; mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3, dt_period=1e-3,
            parallel::Bool=true, distributed::Bool=true, uneven::Bool=false, write_file::Bool=true, filename::Union{String,Int}=-1, optional_conditions::NamedTuple=NamedTuple())
        pairs = initialize_pairs(sim_conditions.M_WDs, sim_conditions.M_NSs, optional_conditions)
        if parallel
            if distributed
                pairs_distr = distribute(pairs)
                @distributed for pair in pairs_distr
                    M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc = pair
                    optional = [mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc]
                    integrate(M_WD; M_NS, mu_e, 
                        orbit_limit, step_limit, mass_limit, dt_period,
                        parallel, uneven, write_file, filename, optional)
                    GC.gc() # Free up memory
                end    
                # Need to be careful here because apparently having each process read and write to disk can worsen performance if not done carefully
            else
                @floop ThreadedEx(basesize = 8) for pair in pairs
                    M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc = pair
                    optional = [mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc]
                    integrate(M_WD; M_NS, mu_e, 
                        orbit_limit, step_limit, mass_limit, dt_period,
                        parallel, uneven, write_file, filename, optional)
                    GC.gc() # Free up memory
                end
            end
        else
            for pair in pairs
                integrate(M_WD; M_NS, mu_e, 
                    orbit_limit, step_limit, mass_limit, dt_period,
                    parallel, uneven, write_file, filename, optional)
                GC.gc() # Free up memory
            end
        end
    end

    """
        Helper function for distributed memory parallelization mode. Makes distributed pairs for all inputted parameter arguments specified by the user, including if the user specified to iterate over the optional parameter arguments.
    """
    function initialize_pairs(M_WDs, M_NSs, optional_conditions)
        mu_es = hasproperty(optional_conditions, :mu_es) ? optional_conditions.mu_es : (optional_default[1],)
        equations = hasproperty(optional_conditions, :WD_mass_radius_equations) ? optional_conditions.WD_mass_radius_equations : (optional_default[2],)
        As = hasproperty(optional_conditions, :As) ? optional_conditions.As : (optional_default[3],)
        poly_shells = hasproperty(optional_conditions, :use_polytropic_density) ? optional_conditions.use_polytropic_density : (optional_default[4],)
        a_ring_fracs = hasproperty(optional_conditions, :a_ring_fracs) ? optional_conditions.a_ring_fracs : (optional_default[5],)
        modes = hasproperty(optional_conditions, :mass_loss_modes) ? optional_conditions.mass_loss_modes : (optional_default[6],)
        eta_accs = hasproperty(optional_conditions, :eta_accs) ? optional_conditions.eta_accs : (optional_default[7],)
        period_calcs = hasproperty(optional_conditions, :use_period_calc) ? optional_conditions.use_period_calc : (optional_default[8],)

        for mass_loss_mode in modes
            first_char = lowercase(first(mass_loss_mode))
            @assert first_char == "j" || first_char == "i" || first_char == "c"
        end

        # Apologies for the syntax below, it isn't elegant but it makes the code run
        pairs = [(M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc) for M_WD in M_WDs for M_NS in M_NSs for mu_e in mu_es for equation in equations for A in As for poly_shell in poly_shells for a_ring_frac in a_ring_fracs for mode in modes for eta_acc in eta_accs for period_calc in period_calcs if M_WD < M_NS]
        return pairs
    end
end