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
    #nu=[a,M_WD,M_NS,R_WD,J,theta,t]
    #theta=[a,M_WD,M_NS,R_WD,J,theta]
    function dnu_dr!(du::Vector{Float64}, nu::Vector{Float64},p::Vector{Any},t::Float64; parallel::Bool=true, uneven::Bool=false, optional=optional_default)
        theta=nu[1:6]
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
    
        dt_dr=(dtheta[1])^(-1)

        dnu_dr=[1, dtheta[2] * dt_dr, dtheta[3] * dt_dr, dtheta[4] * dt_dr, dtheta[5] * dt_dr, dtheta[6] * dt_dr, dt_dr]
        # dt/da, dM_WD/da, dM_NS/da, dR_WD/da, dtheta/da

    
        #println(t_dtheta_dr)

        du.=dnu_dr

    end

    
        

    """ Normal parallelization implementation of in-place acceleration! function""" 
    function dnu_dr_even!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        dnu_dr!(du, theta, p, t; optional=optional)
    end

        
    """ Serial implementation of in-place acceleration! function"""
    function dnu_dr_serial!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        dnu_dr!(du, theta, p, t; parallel=false, optional=optional)
    end


    function integrate_r_0(nu,sign,
                            optional=optional_default;
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9())

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
                        dnu_dr_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dr_serial!(du, θ, p, t; optional=optional)
            end
        end


        if sign==1
            final_r=nu[1]*10
        else
            final_r=0
        end

        
        prob = ODEProblem(t_acceleration!, nu, (nu[1], final_r), Any[])
        integ = init(prob,alg)

        Flag=true

        while Flag
            try 
                step!(integ)
                println(integ.t)
                println(integ.u)
            catch
                Flag=false
            end
        end

        sol=solve(prob,alg,dt=1e-2,adaptive=false,save_everystep=false)
        print(sol.u)
        
    end

    function integrate_r_1(nu,sign;
                            optional=optional_default,
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9(),dr=0.1,
                            frac=5,step_count=5,
                            res=1e-4)

        @assert sign == 1 || sign == -1 #posive if binary is widening, negative if binary is hardening

        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.


        println(dr)
    
        t_acceleration! = nothing
        redirect_stderr(devnull) do
            if uneven
                throw("We don't do that here.") 
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dr_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dr_serial!(du, θ, p, t; optional=optional)
            end
        end


        flag=true
        out=zeros(1000,7)
        i=0

        steps=10*2^step_count
    
        while flag
            try
                prob = ODEProblem(t_acceleration!, nu, (nu[1], nu[1]+sign*steps*dr), Any[])
                sol=solve(prob,alg,dt=dr,adaptive=false,save_everystep=true)
                #println(sol.u)
            
                for j in i*10+1:i*11
                    out[i*10+1+j,:].=sol.u[j]
                end
                nu.=sol.u[end]
                println(nu)
                
            catch
                if steps>10
                    steps/=2
                else
                    flag=false
                    if dr>res
                        new_out=integrate_r_1(nu,sign;
                                optional=optional,
                                 parallel=parallel,
                                 uneven=uneven,
                                 alg=alg,dr=dr/frac)
                        println(1:size(new_out)[1])
                        for j in 1:(size(new_out)[1]) 
                            out[i*10+1+j,:].=new_out[j]
                        end
                        return(out[1:i*10+1+size(new_out),:])
                    else
                        return(out[1:i*10+1,:])
                    end
                end
            end
        end
    end


    """Integrate the change from the initial state rather than the initial state itself to mitigate floating point errors """
    function dnu_dr_sneaky(nu_init::Vector{Float64})
            function ftn!(du::Vector{Float64}, delta_nu::Vector{Float64},p::Vector{Any},t::Float64; parallel::Bool=true, uneven::Bool=false, optional=optional_default)
                 dnu_dr!(du, nu_init.+delta_nu,p,t; parallel=parallel, uneven=uneven, optional=optional)
            end
            return ftn!
        end



    """ Normal parallelization implementation of in-place acceleration! function""" 
    function dnu_dr_even_skeaky!(du::Vector{Float64}, nu_init::Vector{Float64}, delta_nu::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        ftn!=dnu_dr_sneaky(nu_init)
        ftn!(du, delta_nu, p, t; optional=optional)
    end

        
    """ Serial implementation of in-place acceleration! function"""
    function dnu_dr_serial_skeaky!(du::Vector{Float64}, nu_init::Vector{Float64}, delta_nu::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        ftn!=dnu_dr_sneaky(nu_init)
        ftn!(du, delta_nu, p, t; parallel=false, optional=optional)
    end


    function integrate_r_sneaky(nu_init,sign;
                            delta=zeros(7),
                            optional=optional_default,
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9(),dr=0.1,
                            frac=5,step_count=5,
                            res=1e-4)

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
                        dnu_dr_even_skeaky!(du, nu_init, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dr_serial_skeaky!(du, nu_init, θ, p, t; optional=optional)
            end
        end


        flag=true
        out=zeros(1000000,7)
        i=0

        steps=10*2^step_count
        println(dr)
    
        while flag
            try
                prob = ODEProblem(t_acceleration!, delta, (delta[1], delta[1]+sign*steps*dr), Any[])
                sol=solve(prob,alg,dt=dr,adaptive=false,save_everystep=true)
                #println(sol.u)
            
                for j in i*10+1:i*11
                    out[i*10+1+j,:].=sol.u[j]
                end
                delta.=sol.u[end]
                println(delta)
                i+=1
                
            catch
                if steps>10
                    steps/=2
                else
                    flag=false
                    if dr>res
                        new_out=integrate_r_sneaky(nu_init,sign;
                            delta=delta,
                            optional=optional,
                             parallel=parallel,
                             uneven=uneven,
                             alg=alg,dr=dr/frac)
                        println(1:size(new_out)[1])
                        #=for j in 1:(size(new_out)[1]) 
                            out[i*10+1+j,:].=new_out[j]
                        end
                        return(out[1:i*10+1+size(new_out),:])=#
                    else
                        return(out[1:i*10+1,:])
                    end
                end
            end
        end
    end

    







    function make_one_step_dtheta_dr(sign,
                            optional=optional_default;
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9())
    
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
                        dnu_dr_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        dnu_dr_serial!(du, θ, p, t; optional=optional)
            end
        end
    
        # Dummy initial state (will be overwritten by reinit!)
        θ0 = ones(7)
    
        prob = ODEProblem(t_acceleration!, θ0, (0.0, 1.0), Any[])
        println("hi")
        integ = init(prob, alg;
                     save_everystep=false,
                     save_start=false,
                     save_end=false)
        println("Bye")
    
        function ftn!(nu::Vector{Float64};
                      dr_frac::Float64=1e-4,
                      dr::Float64=-99.)::Tuple{Float64,Vector{Float64}}
    
            @assert length(nu) == 7
    
            # Determines separation step as a fraction of the separation
            if dr == -99
                @inbounds begin 
                    dr = sign * dr_frac * nu[1] 
                    #sign insures that dt<0 if the binary is hardening
                    #or dt>0 if the binary is widening
                end
            end
    
            # Reset integrator state 
            reinit!(integ; u0=nu, t0=0.0) 
            #no explicit time-dependent variable here 
            #because true time is inside nu 
            #so it is okay to start over at t0=0
            #(t0 does not mean anything)
    
            # Force integration to exactly dr
            empty!(integ.opts.tstops)
            add_tstop!(integ, dr)
    
            advance_to_tstop!(integ)
    
            nu .= integ.u
            return copy(nu)
        end
    
        return ftn!
    end
        


    function integrate_r(nu, sign; 
                                mu_e::Float64=2.,
                                step_limit=1e4,
                                mass_limit::Float64=1e-2, dr_frac=1e-4,
                                parallel::Bool=true, uneven::Bool=false,
                                write_file::Bool=true, filename::Union{String,Int}=-1,
                                optional=optional_default,
                                checkpoint_every::Int=500,
                                verbose::Bool=false)


        println(nu)
        println("UPDATED")
    
        one_step! = make_one_step_dtheta_dr(sign, optional; parallel=parallel, uneven=uneven) # NEW

        println(nu[1:6])
        ddI=Math.ddI_from_theta(nu[1:6],optional=optional)
        println("ddI",ddI)

        datumses=zeros(step_limit+1,8)
        datumses[1,:].=[sqrt(ddI[1,1]^2+ddI[1,2]^2), nu... ]
        i=1

        while i<step_limit && nu[2]>mass_limit
            println(dr_frac)
            println(one_step!)
            dt,nu=one_step!(nu; dr_frac=dr_frac)
            ddI=Math.ddI_from_theta(nu[1:6]; optional=optional)
            i+=1
            datumses[i,:].=[sqrt(ddI[1,1]^2+ddI[1,2]^2), nu... ]
            if mod(i,10)==0 && verbose
                println(i," ",  theta)
            end
        end


        datumses=datumses[1:i,:]
    
        if write_file #saves datumses to a csv in directory "Output"
            Save.save_as_csv(datumses,M_WD,M_NS,optional.A,optional.mode,filename=filename)
        end
        
        return(datumses)
    end

end