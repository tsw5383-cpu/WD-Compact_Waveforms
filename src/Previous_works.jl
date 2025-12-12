#Here is where I am putting outdated code no longer used in our modules for purposes of documentation.
#(Even though I know that GitHub should save everything, I am still nervous to delete old code, so I am putting it here.)




#------------------------------------------------------------------------------------------------------------------------
#From polytrope.jl:
#=
    """ Returns a function f(xi, theta) where theta = [θ, θ′]
    #Computes Lane–Emden 2nd-order equation in first-order form."""
    @memoize function polytrope_a(n)
        function f!(dθ, θ, p, xi)
            # θ = [theta, dtheta/dxi]
            # dθ = [d theta/dxi, d² theta/dxi²]
            # Lane-Emden: (1/xi^2) d/dxi(xi^2 θ′) = -θ^n
            θ_0=max(θ[1],0) #avoid negative exponentiation
            a = -2 * θ[2] / xi - θ_0^n
            dθ[1] = θ[2]
            dθ[2] = a
        end
        return f!
    end
    
    # ----------------------------
    # lane_emden_int(xi, sol)
    # ----------------------------
    """continuous approximation for dimensionless polytrope density"""
    @memoize function lane_emden_int(xi, sol)
        if xi > sol.t[end]
            return 0.0
        end
        θ = sol(xi)       # sol(x) gives interpolated vector [θ, θ′]
        return θ[1]
    end
    
    """ dm/dxi = xi^2 * θ^n  (standard Lane–Emden mass integrand)"""
    function dm_dxi(xi, sol, n)
        θ = lane_emden_int(xi, sol)
        return xi^2 * θ^n
    end
    
    
    """returns a table for a scaled radius and mass contained within
    the radius"""
    @memoize function scaled_polytrope_table(n)
    
        f! = polytrope_a(n)
    
        # Solve Lane-Emden from xi = 0.001 to 4 if n<=1.5
        #or up to xi=10 elsewise
    
        θ0 = [1.0, 0.0]
        if n<=1.5
            tspan = (0.001, 4.0)
        else
            tspan = (0.001, 10.0)
        end
        prob = ODEProblem(f!, θ0, tspan, nothing)
        sol = solve(prob; reltol=1e-9, abstol=1e-9, dtmax=1e-3)  # RK45 default

        # Compute enclosed mass at 101 radii
        xis = range(0, sol.t[end], length=101)
    
        m_encl = [ quadgk(x -> dm_dxi(x, sol, n), 0, xi)[1] for xi in xis ]
    
        # Normalize radius and mass
        total_mass = quadgk(x -> dm_dxi(x, sol, n), 0, sol.t[end])[1]
    
        scaled_polytrope = zeros(101, 2)
        scaled_polytrope[:, 1] .= xis ./ sol.t[end]     # normalized radius
        scaled_polytrope[:, 2] .= m_encl ./ total_mass  # normalized mass
    
        return scaled_polytrope
    end

    """Reads polytrope solutions from a pre-calculated table so that they do not 
    need to be recalculated each time or writes such a file if it has not yet been 
    written. """
    @memoize function load_polytrope_table(n;rounding=2)
        n=round(n;digits=rounding)
        filename = "Polytrope_Tables_n=$(n).csv"
        try
            return Matrix(CSV.read(filename, DataFrame))
        catch
            tab = scaled_polytrope_table(n)
            CSV.write(filename, DataFrame(tab))
            return tab
        end
    end

    #scaled_poly_15=load_polytrope_table(1.5)
    #scaled_poly_30=load_polytrope_table(3.0)
    

    """ For a white dwarf with radius R_WD, mass M_WD, and polytropic index 
    n, returns the mass outside of a given radius, R_RL"""
    @memoize function mass_outside_radius(R_WD,R_RL,M_WD,n)

        if R_RL>R_WD
            return 0
        end
        
        if n==1.5
            tab=scaled_poly_15
        elseif n==3.0
            tab=scaled_poly_30
        else
            tab=load_polytrope_table(n)
        end

        r=R_RL/R_WD

        return (M_WD*(1-linear_interpolation(tab[:,1], tab[:,2])(r)))
    end

=#

#___________________________________________________________________________________________________________________________________
#math.jl
#=
function time_derivative(ftn::Function, ind_vars_and_derivatives::AbstractVector{T}; optional=optional_default) where {T<:Real}
        #ftn is the function to be taken a time derivative of
        #ind_vars is an array with the values of the arguments of ftn
        #derivatives is an array with the known time derivatives of ind_vars
        @assert mod(length(ind_vars_and_derivatives),6)==0
    
    	ind_vars=@view ind_vars_and_derivatives[1:end-6] 
        #There may be some issues with using view here due to the way that ForwardDiff uses vectors,
        #so I think it is better to use a full vector in this case
        ind_vars= ind_vars_and_derivatives[1:end-6]
    	derivatives=@view ind_vars_and_derivatives[7:length(ind_vars_and_derivatives)]

        #defining t_ftn to be the function with optional parameters fixed to input optional
        function t_ftn(vars)
            return(ftn(vars; optional=optional))
        end
    
    	# ftn_at_ind_vars is the function evaluated at the parameters in ind_vars
        ftn_at_ind_vars = DiffResults.GradientResult(ind_vars)
    	ForwardDiff.gradient!(ftn_at_ind_vars,t_ftn,ind_vars)

        # DiffResults.grad(ftn_at_ind_vars) gives the gradient of the function, evaluated at the parameters in ind_vars
        return(sum(DiffResults.gradient(ftn_at_ind_vars) .* derivatives))#chain rule
    end=#

#-----------------------------------------------------------------------------------------------------------------------------
#Integrate.jl



#OLD CODE
    #PROBLEM: This function creates a new ODE problem and solver every time and is thus terribly memory inefficient.
    
    #=function one_step!(theta::Vector{Float64}; optional=optional_default, dt_period::Float64=1e-2, dt::Float64=-1., parallel::Bool=true, uneven::Bool=false)::Tuple{Float64, Vector{Float64}}
        @assert length(theta)==6
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

        t_acceleration! = nothing
    
        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        redirect_stderr(devnull) do 
            if uneven
                prob = ODEProblem(Math.acceleration_uneven!,theta,(0.,dt),[])
            else
                t_acceleration! = parallel ? 
                    (du, θ, p, t) -> Math.acceleration_even!(du, θ, p, t; optional=optional) :
                    (du, θ, p, t) -> Math.acceleration_serial!(du, θ, p, t; optional=optional) 
                prob = ODEProblem(t_acceleration!, theta, (0., dt), []) #THIS IS THE PROBLEM HERE
            end
        end
        theta .= solve(prob, alg=Vern9()).u[end]
        dt,theta
    end=#



#REPLACEMENT FOR OLD one_step

    #=#hyperparameters which will be used to define the acceleration function
    mutable struct Hyperparams
        optional::Vector{Any}
        parallel::Bool
        uneven::Bool
    end


    """ """
    function t_accel!(du::Vector{Float64}, θ::Vector{Float64}, p::Hyperparams, t::Float64, u::Vector{Any})
        if p.uneven
            Math.acceleration_uneven!(du, θ, u, t)
        elseif p.parallel
            Math.acceleration_even!(du, θ, u, t; optional=p.optional)
        else
            Math.acceleration_serial!(du, θ, u, t; optional=p.optional)
        end
    end
    #=function t_accel!(du, θ, p::Hyperparams, t)
        if p.uneven
            Math.acceleration_uneven!(du, θ, nothing, t)
        elseif p.parallel
            Math.acceleration_even!(du, θ, nothing, t; optional=p.optional)
        else
            Math.acceleration_serial!(du, θ, nothing, t; optional=p.optional)
        end
    end
        #=else
           t_acceleration! = p.parallel ? 
                    (du, θ, u, t) -> Math.acceleration_even!(du, θ, u, t; optional=optional) :
                    (du, θ, u, t) -> Math.acceleration_serial!(du, θ, u, t; optional=optional)
            t_acceleration!
        end
    #end=#=#
    
    """
    Creates and returns a function `step!(theta, dt)` that advances `theta` by dt and returns dt,theta.
    This reuses an ODEIntegrator internally for speed and to avoid global captures as in our old approach.
    """
    #=function make_one_step(optional; parallel::Bool=true, uneven::Bool=false, alg=Vern9())
        p = Hyperparams(optional, parallel, uneven)
    
        # build a template problem with a *dummy* initial condition
        u0 = zeros(6)                
        prob = ODEProblem(t_accel!, u0, (0.0, 1.0), p)  
    
        # create an integrator object that we will reuse
        integ = init(prob, alg; save_everystep=false)
    
        # returned function: advance theta by dt
        function t_step!(theta::Vector{Float64}; dt_period::Float64=1e-2, dt_in::Float64=-1.0)
            @assert length(theta)==6
    
            dt = dt_in
            if dt == -1.0
                @inbounds begin
                    a = theta[1]
                    M_WD = theta[2]
                    M_NS = theta[3]
                    dt = dt_period * Math.period(a, M_WD + M_NS)
                end
            end
    
            # redefine initial conditions to our problem
            integ.u .= theta
    
            # advances integrator by dt
            step!(integ, dt)   
    
            # copy result back to theta (in-place)
            theta .= integ.u
    
            return dt, theta
        end
    
        return t_step!
    end=#

    function make_one_step(optional::Vector{Any};
        parallel::Bool=true,
        uneven::Bool=false,
        alg=Vern9(),
        tspan=(0.0, 1.0)
    )
    
        # Build parameter container
        p = Hyperparams(parallel, uneven, optional)
    
        # ODE RHS
        f = (du, θ, t) -> t_accel!(du, θ, p, t)
    
        # Build a problem with dummy initial condition
        prob = ODEProblem(f, zeros(9), tspan, p)
    
        # Build integrator
        integrator = init(prob, alg; dt=1e-2)
    
        # Return stepping function
        return function(theta; dt_period=1e-2, dt_in=-1)
            dt = dt_in
            if dt == -1.0
                @inbounds begin
                    a = theta[1]
                    M_WD = theta[2]
                    M_NS = theta[3]
                    dt = dt_period * Math.period(a, M_WD + M_NS)
                end
            end
        
            # overwrite integrator state with user-supplied theta
            integrator.u .= theta
    
            # perform one step of size dt_period
            step!(integrator, dt)
    
            # return updated dt and theta
            return dt, copy(integrator.u)
        end
    end=#