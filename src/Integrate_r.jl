module Integrate_r
    using Distributed
    using OrdinaryDiffEq
    using FLoops
    using DistributedArrays
    include("math.jl")
    include("save.jl")
    include("Integrate.jl")

    const optional_default = (mu_e = 2.0, equation = 1, A = 10.0, poly_shell = false, a_ring_frac = 1.0, mode = :isotropic, eta_acc = 0.1, period_calc = true)

    #calculates dtheta/dr
    function dtheta_dr!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; parallel::Bool=true, uneven::Bool=false, optional=optional_default)
        if parallel && uneven
            #keyword argument dictates the style of parallelization to be used 
            #(mostly for benchmarking purposes)
            dtheta,t_ddI,t_dddI=Math.ddI_dddI_from_theta_uneven(theta)
        else
            dtheta,t_ddI,t_dddI=Math.ddI_dddI_from_theta(theta; parallel=parallel, optional=optional)
        end
    
        #forces due to RR
        t_dJ_rr::Float64=Math.dJ_rr(t_ddI,t_dddI)
        @inbounds begin
            a::Float64 = theta[1]
            J::Float64 = theta[5]
            t_da_rr::Float64 = 2.0 * a * t_dJ_rr / J

            #println(t_ddI,t_dddI)
            #println(t_dJ_rr)
            #println(t_da_rr)# For debugging--DELETE LATER
            #println()
        
            #adding RR forces to MT forces
            dtheta[1]+=t_da_rr
            dtheta[5]+=t_dJ_rr
        end


        println(dtheta)
    
        dt_dr=(dtheta[1])^(-1)

        t_dtheta_dr=[dt_dr, dtheta[2] * dt_dr, dtheta[3] * dt_dr, dtheta[4] * dt_dr, dtheta[5] * dt_dr, dtheta[6] * dt_dr]

        #println(t_dtheta_dr)

        du.=t_dtheta_dr

    end

    """ Normal parallelization implementation of in-place acceleration! function""" 
    function dtheta_dr_even!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        dtheta_dr!(du, theta, p, t; optional=optional)
    end

        
    """ Serial implementation of in-place acceleration! function"""
    function dtheta_dr_serial!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        dtheta_dr!(du, theta, p, t; parallel=false, optional=optional)
    end

    function make_one_step_dtheta_dr(optional=optional_default;
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9())

        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        t_acceleration! = nothing
        redirect_stderr(devnull) do
            if uneven
                raise 
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) ->
                        dtheta_dr_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        dtheta_dr_serial!(du, θ, p, t; optional=optional)
            end
        end
    
        # Dummy initial state (will be overwritten by reinit!)
        θ0 = zeros(6)
    
        prob = ODEProblem(t_acceleration!, θ0, (0.0, 1.0), Any[])
        integ = init(prob, alg;
                     save_everystep=false,
                     save_start=false,
                     save_end=false)
    
        function ftn!(theta::Vector{Float64};
                      dr_frac::Float64=1e-4,
                      dt::Float64=-1.0)::Tuple{Float64,Vector{Float64}}
    
            @assert length(theta) == 6
    
            # Determines separation step as a fraction of the separation
            if dr == -1.0
                @inbounds begin 
                    dt = dt_period * theta[1]
                end
            end
    
            # Reset integrator state 
            reinit!(integ; u0=theta, t0=0.0) 
            #no explicit time-dependent variable here 
            #(i.e., time translation invariant) 
            #so it is okay to start over at t=0 
    
            # Force integration to exactly dt
            empty!(integ.opts.tstops)
            add_tstop!(integ, dt)
    
            advance_to_tstop!(integ)
    
            theta .= integ.u
            return dt, copy(theta)
        end
    
        return ftn!
    end
        



end