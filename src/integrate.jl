module Integrate
include("./core_math.jl")

"""Defines theta at initial Roche Lobe overflow for  """
function initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2., mode::String="Debug")::Vector{Float64}
    #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
    #M_NS=1.35 M_\odot is typical NS mass
    if mode=="Debug"
        @assert M_WD<1.4 #Chandrasekar mass, beyond which we cannot calculate R_WD(see initial function)
    end
    a=RL_contact(M_WD,M_WD,mu_e=mu_e)
    J=circular_J(a,M_WD,M_WD)
    R_WD=R0_WD(M_WD)
    return([a, M_WD, M_NS, R_WD, J, 0.])
end

#"Below is the heart of the integration process"
"""Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken, the updated theta parameters, and the instantaneous second derivative of the quadupole moment at the end of the timestep."""
function one_step(theta::Vector{Float64}; dt_period::Float64=1e-2, dt::Float64=-1., mode::String="Debug")::Tuple{Float64, Vector{Float64}, Array{Float64,2}}
    if mode=="Debug"
        @assert length(theta)==6
    end
    #integrates forward one step with time interval either given explicitly (dt!=-1)
    #or as a fraction of the orbital period determined by dt_period
    if dt==-1
        dt=dt_period*period(theta[2]+theta[3], theta[1])
        #print(dt)
    end
    prob = ODEProblem(acceleration!,theta,(0,dt),[])
    updated_theta = solve(prob, alg=Vern9()).u[end]

    #=Calculating ddI at the end of the step=#
        
    t_ddI=ddI_from_theta(updated_theta)
    #println(updated_theta)
    return(dt,updated_theta,t_ddI)
end

"""Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. The function also saves the run to a file named by the date and initial conditions if `write_file`. """
function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
        orbit_limit::Float64=1e5, step_limit::Int=10^7,
               write_file::Bool=true, filename::Union{String,Int}=-1)::AbstractMatrix{<:Real}
    #returns time series second derivatives of the quadrupole moment
    t=0.
    theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
    ddI=ddI_from_theta(theta)
    #ts=zeros(step_limit)
    #ddIs=zeros(step_limit,2,2)
    datumses=zeros(step_limit+1,9)
    datumses[1,:].=[t, ddI[1,1], ddI[1,2], theta... ]
    #print(ts_plus_times)
    i=1
    
    Flag=true
    
    while theta[6]<2*π*orbit_limit && i<step_limit && theta[2]>0 && Flag
        #integrate until time or step limit are elapsed or WD completely dissipates
        #println(theta)
        #println(isa(theta,Vector{Float64}))
        try
            dt,t_theta,t_ddI=one_step(theta; dt_period=1e-2)
            i+=1
            t+=dt
        
            datumses[i,:].=[t, t_ddI[1,1], t_ddI[1,2], t_theta... ]
    
            theta=t_theta
        catch
            Flag=false
        end
    end

    datumses=datumses[1:i,:]

    if write_file #saves datumses to a csv in directory "Output"
        save_as_csv(datumses,M_WD,M_NS,filename=filename)
    end
    
    return(datumses)
end