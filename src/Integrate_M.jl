module Integrate_M
    using Distributed
    using OrdinaryDiffEq
    using FLoops
    using DistributedArrays
    using Printf
    include("math.jl")
    include("save.jl")
    include("Integrate.jl")

    const optional_default = (mu_e = 2.0, equation = 1, A = 10.0, poly_shell = false, a_ring_frac = 1.0, mode = :isotropic, eta_acc = 0.1, period_calc = true)

    #calculates dtheta/dM_WD
    #nu=[a,M_WD,M_NS,R_WD,J,theta,t]
    #theta=[a,M_WD,M_NS,R_WD,J,theta]
    function dnu_dM!(du::Vector{Float64}, nu::Vector{Float64},p::Vector{Any},t::Float64; parallel::Bool=true, uneven::Bool=false, optional=optional_default,evolve_R_WD=true)
        theta=nu[1:6]
        if !evolve_R_WD #if t_dyn<<time step, the radius of the white dwarf evolves essentially instantanously 
            theta[4]=Math.R0_WD(theta[2])
        end
        #println("theta",theta)
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


        #println(dtheta)
    
        dt_dr=(dtheta[2])^(-1) #Actually dt_dM_WD (just easier to call it this because copy-and-paste)

        if evolve_R_WD #evolve R_WD on dynamical timescale
            dnu_dr=[dtheta[1] * dt_dr, 1, dtheta[3] * dt_dr, dtheta[4] * dt_dr, dtheta[5] * dt_dr, dtheta[6] * dt_dr, dt_dr]
        else #assume R_WD is R_WD_0 (i.e., equilibrium radius)
            dnu_dr=[dtheta[1] * dt_dr, 1, dtheta[3] * dt_dr, 0, dtheta[5] * dt_dr, dtheta[6] * dt_dr, dt_dr]
        end
        # dt/da, dM_WD/da, dM_NS/da, dR_WD/da, dtheta/da

    
        #println(t_dtheta_dr)

        du.=dnu_dr

    end

    
        

    """ Normal parallelization implementation of in-place acceleration! function""" 
    function dnu_dM_even!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default,evolve_R_WD=false)
        dnu_dM!(du, theta, p, t; optional=optional,evolve_R_WD=evolve_R_WD)
    end

        
    """ Serial implementation of in-place acceleration! function"""
    function dnu_dM_serial!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default,evolve_R_WD=false)
        dnu_dM!(du, theta, p, t; parallel=false, optional=optional,evolve_R_WD=evolve_R_WD)
    end





    """Integrate the change from the initial state rather than the initial state itself to mitigate floating point errors """
    function dnu_dM_sneaky(nu_init::Vector{Float64})
            function ftn!(du::Vector{Float64}, delta_nu::Vector{Float64},p::Vector{Any},t::Float64; parallel::Bool=true, uneven::Bool=false, optional=optional_default,evolve_R_WD=false)
                 dnu_dM!(du, nu_init.+delta_nu,p,t; parallel=parallel, uneven=uneven, optional=optional,evolve_R_WD=evolve_R_WD)
            end
            return ftn!
        end



    """ Normal parallelization implementation of in-place acceleration! function""" 
    function dnu_dM_even_sneaky!(du::Vector{Float64}, nu_init::Vector{Float64}, delta_nu::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default,evolve_R_WD=evolve_R_WD)
        (dnu_dr_sneaky(nu_init))(du, delta_nu, p, t; optional=optional,evolve_R_WD=evolve_R_WD)
    end

        
    """ Serial implementation of in-place acceleration! function"""
    function dnu_dM_serial_sneaky!(du::Vector{Float64}, nu_init::Vector{Float64}, delta_nu::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default,evolve_R_WD=false)
        (dnu_dr_sneaky(nu_init))(du, delta_nu, p, t; parallel=false, optional=optional,evolve_R_WD=evolve_R_WD)
    end



    """Integrates forward to a specific mass (M) in the M_WD-basis"""
    function integrate_M_specific(nu_init,sign; M=0.05
                            delta=zeros(7),
                            optional=optional_default,
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9(),dr=0.1,
                            frac=5,step_count=8,
                            res=1e-4,evolve_R_WD=true,
                            final_return=true)

        @assert sign == 1 || sign == -1 #posive if binary is widening, negative if binary is hardening

        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        t_acceleration! = nothing
        redirect_stderr(devnull) do
            if uneven
                throw("We don't do that here.") 
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dM_even_sneaky!(du, nu_init, θ, p, t; optional=optional,evolve_R_WD=evolve_R_WD) :
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dM_serial_smeaky!(du, nu_init, θ, p, t; optional=optional,evolve_R_WD=evolve_R_WD)
            end
        end


        prob = ODEProblem(t_acceleration!, delta, (0, M-nu_init[2]), Any[]) 
                
        sol=solve(prob,alg,dt=dr,alg,reltol=1e-9,adaptive=true,save_everystep=true)

        return(sol)
    
    end

end