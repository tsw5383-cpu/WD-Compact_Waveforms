begin
	using LinearAlgebra
	using BenchmarkTools
	using Profile, ProfileCanvas 
	using PlutoUI, PlutoTeachingTools, PlutoLinks
	using Test

    using Pkg
    function ensure_package_installed(pkg_name::String)
        # Check if the package is installed
        installed_packages = Pkg.dependencies()
        if !any(p -> p.name == pkg_name, values(installed_packages))
            println("Package '$package_name' not found. Installing...")
            Pkg.add(package_name)
            println("Package '$package_name' installed successfully.")
        end
    end
    ensure_package_installed("OrdinaryDiffEq")
    ensure_package_installed("ForwardDiff")
    ensure_package_installed("SciMLBase")
    ensure_package_installed("CSV")
    ensure_package_installed("DataFrames")
    ensure_package_installed("Dates")
    ensure_package_installed("DelimitedFiles")
    ensure_package_installed("FilePathsBase")

	using ForwardDiff, SciMLBase
    using OrdinaryDiffEq
	#packages used for chain rules of large derivatives
	#and integrating differential equations

	using CSV, DataFrames
	using Dates,DelimitedFiles, Printf
	using FilePathsBase: joinpath
	#for saving data as files
	
	current_path=@__DIR__
	include(current_path*"/victor-ftns_only_no_types.jl")

	const G::Float64 = 1.327e11   # km^3 / M⊙ / s^2
	const c::Float64 = 2.998e5    # km/s
	#victor-functions_only.jl:
	#file with a lot of helpful functions
	#This is the biggest current issue: 
	#I do not know how to use the functions 
	#in victor.jl in tristan.jl. This causes 
	#us to fail our final integration test.
	#EDIT: I figured this out after much 
	#bashing of my head against the wall.
	#EDIT: Removing types for compatibility
	#with ForwardDiff.
	
end;

# *** #
"""Smooth approximation to the minimum function"""
function Boltzmann_Max(x1::T, x2::T; α=-10.) where T<:Real
    e_α_x1 = exp(α*x1)
    e_α_x2 = exp(α*x2)
    return((x1*e_α_x1+x2*e_α_x2)/(e_α_x1+e_α_x2))
end

begin
	""" Smooth accretion rate of the neutron star """
	function dM_NS(a::T, R_WD::T, M_WD::T, M_NS::T;
			A=10., mu_e=2.0, eta_acc=0.1, α=-10.) where T<:Real
	    t_dM_WD=decretion_rate(a,a,R_WD,M_WD,M_NS,A=A)
		
	    t_edd=eddington_rate(M_NS, mu_e=mu_e, eta_acc=eta_acc)
		return(min(-t_dM_WD,t_edd))
	end

    # *** #
	"""Reduced mass in effective one body formalism"""
	function mu(M_WD::T, M_NS::T) where T<:Real
		return((M_WD*M_NS)/(M_WD+M_NS))
	end

	# *** #
	""" Time derivative of angular parameter  θ"""
	function dθ(J::T,a::T,M_WD::T, M_NS::T; mode::String="Debug") where T<:Real
        if mode=="Debug"
    		@assert M_WD > 0 # In production-ready code, consider commenting this out
    		@assert M_NS > 0 # In production-ready code, consider commenting this out
    		@assert a > 0 # In production-ready code, consider commenting this out
        end
		return(J/(a^2*mu(M_WD,M_NS)))
	end   
end;

begin
    # *** #
	"""Given polar orbital parameters, calculates I_{xx} """
	function I_xx(M_WD::T, M_NS::T, r::T, θ::T; mode::String="Debug") where T<:Real
        if mode=="Debug"
    	    @assert M_WD > 0 # In production-ready code, consider commenting this out
    		@assert M_NS > 0 # In production-ready code, consider commenting this out
    	    @assert r > 0 # In production-ready code, consider commenting this out
        end
	    return (mu(M_WD, M_NS) * (x(r,θ))^2)
	end           

    # *** #
	"""Given polar orbital parameters, calculates I_{xy} """
	function I_xy(M_WD::T, M_NS::T, r::T, θ::T; mode::String="Debug") where T<:Real
        if mode=="Debug"
    	    @assert M_WD > 0 # In production-ready code, consider commenting this out
    		@assert M_NS > 0 # In production-ready code, consider commenting this out
    	    @assert r > 0 # In production-ready code, consider commenting this out
        end
	    return mu(M_WD, M_NS) * x(r,θ) * y(r,θ)
	end                
end;

"""Fixes the higher derivatives of the function so that it is formatted in a way that ForwardDiff is comfortable with."""
function fix_higher_derivatives(ftn::Function,higher_derivatives::AbstractVector{T}) where {T} 
	function inner_ftn(gen_theta::AbstractVector{T}) where {T}
		return(ftn(gen_theta,higher_derivatives))
	end
	return ftn 
end;

# *** #
"""Given a function, its vectorized parameters, and the time derivatives of those parameters, returns the time derivative of the function as calculated via the chain rule."""
"""Wraps the vectorized parameters and their derivatives into one input vector ind_vars_and_derivatives"""
function time_derivative(ftn::Function, ind_vars_and_derivatives::AbstractVector{T}; mode::String="Debug") where {T}
    #ftn is the function to be taken a time derivative of
    #ind_vars is an array with the values of the arguments of ftn
    #derivatives is an array with the known time derivatives of ind_vars
	#print(length(ind_vars_and_derivatives))
    if mode=="Debug"
        @assert mod(length(ind_vars_and_derivatives),6)==0
    end

	ind_vars=ind_vars_and_derivatives[1:end-6]
	derivatives=ind_vars_and_derivatives[7:length(ind_vars_and_derivatives)]

	grad=ForwardDiff.gradient(ftn,ind_vars)
    
    return(sum(grad .* derivatives))#chain rule
end

# *** #
"""Given a function, defines its next derivative based on the function's parameters and their time derivatives. `n` determines how many time derivatives of \vec{θ} (i.e., the parameters to be integrated) are parameters of the function itself. """
function next_derivative(in_ftn::Function; n::Integer=1; mode::String="Debug")
    #returns a function which is the next higher derivative of the in_ftn
    #n is the total derivatives past the given analytic form
    if mode=="Debug"
    	@assert n>0
    end
    function ftn(gen_theta_and_derivatives::AbstractVector{T}) where {T}
		#print(in_ftn)
        if mode=="Debug"
            @assert length(gen_theta_and_derivatives)==6*(n+1) 
        end
		#outdated code
		#=gen_theta=gen_theta_and_derivatives[1:end-6]
		gen_dtheta=gen_theta_and_derivatives[end-5:end]=#
        #requires gen_theta to be [theta, dtheta, d...dtheta] for n=number of ds
        #requires gen_theta to be [dtheta, ddtheta d...ddtheta] for n+1=number of ds
        return(time_derivative(in_ftn, gen_theta_and_derivatives))
    end
    return ftn
end

md"The following functions wrap the original derivatives in `victor.jl` in vectors to make life easier when taking derivatives due to mass transfer only. In particular, `theta` is of the form `[a, M_WD, M_NS, R_WD, J, θ]`."
begin
	#=The following functions wrap the original derivatives in vectors to make life easier when taking derivatives due to mass transfer only=#    
	function da_vector(theta::AbstractVector{T}; 
			A=10.0, mu_e=2.0, eta_acc=0.1, 
			a_ring_frac=1.0, mode::String="Isotropic") where {T}

		# promote the numeric keyword args to match the type of the inputs
	    A = T(A)
	    mu_e = T(mu_e)
	    eta_acc = T(eta_acc)
	    a_ring_frac = T(a_ring_frac)

		#unpacking theta
	    a, M_WD, M_NS, R_WD, J, θ=theta

		
	    return(a_dot(a,a,R_WD, M_WD, M_NS;
	            A=A,mu_e=mu_e,eta_acc=eta_acc,
	            a_ring_frac=a_ring_frac,mode=mode))
	end
	
	function dM_WD_vector(theta::AbstractVector{T}; A=10.) where {T}
		A = T(A)
		
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(decretion_rate(a, a, R_WD, M_WD, M_NS, A=A))
	end
	
	function dM_NS_vector(theta::AbstractVector{T};
	        A=10, mu_e=2.0, eta_acc=0.1, α=-10.) where {T}

		A = T(A)
	    mu_e = T(mu_e)
	    eta_acc = T(eta_acc)
	    α = T(α)
		
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(dM_NS(a, R_WD, M_WD, M_NS;
	        A=a, mu_e=mu_e, eta_acc=eta_acc, α=α))
	end
	
	function dR_WD_vector(theta::AbstractVector{T}; mu_e=2.) where {T} 
		mu_e = T(mu_e)
		
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return( dot_R_WD(M_WD, R_WD))
	end
	    
	function dJ_vector(theta::AbstractVector{T}; 
	               A=10.0, mu_e=2.0, eta_acc=0.1, 
	               a_ring_frac=1.0, mode::String="Isotropic") where {T}
		A = T(A)
	    mu_e = T(mu_e)
	    eta_acc = T(eta_acc)
	    a_ring_frac = T(a_ring_frac)
		
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(J_dot(J, a, a, R_WD, M_WD, M_NS, 
	               A=A, mu_e=mu_e, eta_acc=eta_acc, 
	               a_ring_frac=a_ring_frac, mode=mode))
	end
	
	function dθ_vector(theta::AbstractVector{T}) where {T} 
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(dθ(J,a,M_WD, M_NS))
	end
	
	
	#=Same treatment, but for quadrupole moments=#
	function Ixx_vector(theta::AbstractVector{T}) where {T} 
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(I_xx(M_WD, M_NS, a, θ))
	end
	
	function Ixy_vector(theta::AbstractVector{T})  where {T} 
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(I_xy(M_WD, M_NS, a, θ))
	end
end

md"The below block of code uses the vectorized derivatives above to output functions which take higher time derivatives of $\vec{\theta}$ and $I_{ij}$ as a function of $\theta$. "
begin
    dda_vector=next_derivative(da_vector)
    ddM_WD_vector=next_derivative(dM_WD_vector)
    ddM_NS_vector=next_derivative(dM_NS_vector)
    ddR_WD_vector=next_derivative(dR_WD_vector)
    ddJ_vector=next_derivative(dJ_vector)
    ddθ_vector=next_derivative(dθ_vector)

    ddda_vector=next_derivative(dda_vector,n=2)
    dddM_WD_vector=next_derivative(ddM_WD_vector,n=2)
    dddM_NS_vector=next_derivative(ddM_NS_vector,n=2)
    dddR_WD_vector=next_derivative(ddR_WD_vector,n=2)
    dddJ_vector=next_derivative(ddJ_vector,n=2)
    dddθ_vector=next_derivative(ddθ_vector,n=2)

    dI_xx_vector=next_derivative(Ixx_vector)
    ddI_xx_vector=next_derivative(Ixx_vector,n=2)
    dddI_xx_vector=next_derivative(Ixx_vector,n=3)
    
    dI_xy_vector=next_derivative(Ixy_vector)
    ddI_xy_vector=next_derivative(Ixy_vector,n=2)
    dddI_xy_vector=next_derivative(Ixy_vector,n=3)
end;  

#=bundling derivatives together which will prove convenient later=#
begin
    first_derivatives=[da_vector, dM_WD_vector, dM_NS_vector, dR_WD_vector, dJ_vector, dθ_vector]
    second_derivatives=[dda_vector, ddM_WD_vector, ddM_NS_vector, ddR_WD_vector, ddJ_vector, ddθ_vector]
    third_derivatives=[ddda_vector, dddM_WD_vector, dddM_NS_vector, dddR_WD_vector, dddJ_vector, dddθ_vector]
end;

#=Constructs a symmetric 2x2 tensor from a plus and times polarization=#
function tensor_plus_times(plus::Float64, times::Float64)::Matrix{Float64}
    return([plus times; times -plus        
    ])
end

md"The below block of code is central to the interface between the mass transfer and radiation reaction effects. It also sits in the innermost loop of our integrator. This is a good place to peer review."
#=Central Engine of the integrator=#
begin
    # *** #
	"""Calculates the second derivative of the quadrupole moment `t_ddI`  due to kinematics of keplerian motion and mass transfer `dtheta`. Also returns higher derivatives of theta useful for calculating the third derivative of the quadrupole moment."""
	function ddI_helper(theta::AbstractVector{T}) where {T}
	    #calulates a list of values which will be helpful in the next two functions
	    dtheta=zeros(6)
	    #calculating first time derivatives due to MT only
	    @inbounds for i in 1:6
	        dtheta[i]=first_derivatives[i](theta)
	    end

		theta_dtheta=vcat(theta,dtheta)
	    ddtheta=zeros(6)
	    #calculating second time derivatives due to MT only
	    @inbounds for i in 1:6
	        ddtheta[i]=second_derivatives[i](theta_dtheta)
	    end
	
	    theta_ddtheta=vcat(theta_dtheta,ddtheta)
	
	
	    
	    #using orbital elements and first and second time derivatives 
	    #to calculate second derivative of quadrupole moment
	    #(for use in calculating RR acceleration)
	    
	    
	    t_ddI=tensor_plus_times(ddI_xx_vector(theta_ddtheta),
	        ddI_xy_vector(theta_ddtheta)
	    )
	    return(dtheta, ddtheta, theta_ddtheta, t_ddI)
	
	end

	"""Returns the second derivative of the quadrupole moment as calculated in `ddI_Helper`. """
	function ddI_from_theta(theta::Vector{Float64})::Array{Float64,2}
	    return(ddI_helper(theta)[end])
	end
	    
	
	"""Calculates the third derivative of the quadrupole moment and returns it along with the time derivative of theta and the second derivative of the quadrupole moment as calculated in `ddI_Helper`."""
	function ddI_dddI_from_theta(theta::Vector{Float64};ddI_only::Bool=false)::Tuple{Vector{Float64},Array{Float64,2},Array{Float64,2}}
	    #returns time derivatives of theta and second and third derivatives of 
	    #the quadrupole moment due to MT effects
	    
	    dtheta, ddtheta, theta_ddtheta, t_ddI = ddI_helper(theta)
	    #using orbital elements and first through third time derivatives  
	    #to calculate third derivative of quadrupole moment 
	    #(for use in calculating RR acceleration)

		#print(theta_ddtheta)

		dddtheta=zeros(6)
	    #calculating second time derivatives due to MT only
	    @inbounds for i in 1:6
	        dddtheta[i]=third_derivatives[i](theta_ddtheta)
	    end

	    theta_dddtheta=vcat(theta_ddtheta,dddtheta)
	
	    t_dddI=tensor_plus_times(dddI_xx_vector(theta_dddtheta),
	        dddI_xy_vector(theta_dddtheta)
	    )
	    return(dtheta,t_ddI,t_dddI)
	end
end;

# *** #
""" An "in place" implementation of the full time derivative of theta including effects from both mass transfer and radiation reaction.""" 
function acceleration!(du::Vector{Float64}, theta::Vector{Float64},p,t)
    dtheta,t_ddI,t_dddI=ddI_dddI_from_theta(theta)

    #forces due to RR
    t_dJ_rr=dJ_rr(t_ddI,t_dddI)
    @inbounds begin
        t_da_rr=2.0 * theta[1] * t_dJ_rr / theta[5]
    
        #adding RR forces to MT forces
        dtheta[1]+=t_da_rr
        dtheta[5]+=t_da_rr
    end

    #mutate du
    du .= dtheta
	#println(du)
    #return(dtheta)
end;

md"The next several blocks define functions intended to initiate the orbit beginning at mass transfer for initial masses, take an integration step forward calculating waveform, and integrate until the White Dwarf dissipates completely"
begin
	"""Calculates in km the orbital separation for the white dwarf to overflow its Roche Lobe and begin mass transfer"""
	function RL_contact(M_WD::Float64,M_NS::Float64;mu_e::Float64=2.)::Float64
	    #returns separation at which a donor WD overflows its
	    #Roche Lobe
	    R_WD=R0_WD(M_WD,mu_e=mu_e)
	    RL_a1=Roche_Limit(1.,M_WD,M_NS) #roche limit at separation of 1 km
	    return(R_WD / RL_a1)#initial separation in km
	end

	"""Calculates the angular momentum of the system at initial Roche Lobe overflow """	
	function circular_J(a::Float64,M_WD::Float64,M_NS::Float64)::Float64
	    #returns angular momentum for a circular orbit
	    return(( a * G * M_WD ^ 2 * M_NS ^ 2 / (M_WD + M_NS))^(1/2))
	end
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
end

md"Below is the heart of the integration process"
begin
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

	
	"""formats the filename as YYYYMMDD_M_WD=*M_WD*_M_NS=*M_NS*.csv
		if a filename is not provided"""
	function format_filename(M_WD::Float64,M_NS::Float64;filename::Union{String,Int}=-1)::String
		#
		if filename==-1 
		        today_str = Dates.format(Dates.now(), "yyyymmdd")
		        base_filename = "Outputs/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS).csv"
		        filename = base_filename
		        suffix = 2
		        while isfile(filename)
		            filename = "Outputs/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS)_$(suffix).csv"
		            suffix += 1
				end
		end
		return filename
	end
		
	"""writes array of format given in integrate to a file"""
	function save_as_csv(datumses::AbstractMatrix{<:Real},
					  M_WD::Float64,M_NS::Float64; filename::Union{String,Int}=-1)::Nothing
		
		
		filename=format_filename(M_WD,M_NS;filename=filename)
		# Writes the file to a CSV titled filename

		isdir("Outputs") || mkdir("Outputs")
		#creates output directory if not already extant
		
			CSV.write(filename, 
						(t = datumses[:,1],
     					ddI_xx = datumses[:,2],
     					ddI_xy = datumses[:,3],
						a = datumses[:,4],
     					M_WD = datumses[:,5],
     					M_NS = datumses[:,6],
						R_WD = datumses[:,7],
     					J = datumses[:,8],
     					θ = datumses[:,9]))
		return nothing
	end
		
		
		


	"""Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. The function also saves the run to a file named by the date and initial conditions if `write_file`. """
	function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
	        orbit_limit::Float64=1e5, step_limit::Int=10^7,
                   write_file::Bool=true, filename::Union{String,Int}=-1)::AbstractMatrix{<:Real}
	    #returns time series second derivatives of the quadrupole moment
	    t=0.
	    theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
	    #ts=zeros(step_limit)
	    #ddIs=zeros(step_limit,2,2)
		datumses=zeros(step_limit,9)
		#print(ts_plus_times)
	    i=0
		
		Flag=true
	    
	    while theta[6]<2*π*orbit_limit && i<step_limit && theta[2]>0 && Flag
	        #integrate until time or step limit are elapsed or WD completely dissipates
			#println(theta)
			#println(isa(theta,Vector{Float64}))
			try
		        dt,t_theta,t_ddI=one_step(theta; dt_period=1e-2)
		        i+=1
		        t+=dt
				#println(t," ", dt, " ", t_theta[6])
		        
	
				t_ddI=ddI_from_theta(theta)
				#print(t_ddI)
	
				#print([t, t_ddI[1,1], t_ddI[1,2]])
	
				#print(ts_plus_times[i])

				datumses[i,:].=[t, t_ddI[1,1], t_ddI[1,2], t_theta... ]
		
		        theta=t_theta
		        #ddIs[i, :, :].=ddI_from_theta(theta)
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
end;

	
md"Note: For peer review, we do NOT recommend running `integrate` fully (until total disruption) as this may take upward of 10 hours!"

# *** #
"""Given a luminosity distance `DL` in units of Mpc and an orientation parameterized by an inclination `i` and azimuthal angle `φ`, returns a function which calculates the instantaneous gravitational strain waveform based on the quadrupole formula."""
function gravitational_waveform_ftn(i::Float64, φ::Float64, DL::Float64; mode::String="Debug")::Function
    #fixes a function to return gravitational wave given orientation and distance
    #DL is in Mpc
    if mode=="Debug":
        @assert 0<=i<=π
        @assert 0<=φ<2π
        @assert DL<10000 #Ensure correct units for DL
    end
    cosi = cos(i)
    cosφ = cos(φ)
    sini = sin(i)
    sinφ = sin(φ)
    e1=[cosi*cosφ, sini*sinφ, cosi]
    e2=[-sinφ,  cosφ,  0]
    es=[e1,e2]

    matrices=zeros(4,2,2)
    @inbounds for i in 1:2
        @inbounds for j in 1:2
            matrices[2*(i-1)+j ,:,:] = es[i] * es[j] #OUTER PRODUCT OF VECTORS
        end
    end

    factor= 2*G*3.2e-20/c^4 #for conversion from Mpc to km
    
    function ftn(ddI::Vector{Float64})::Vector{Float64}
        #given fixed outer products, calculate strain
        prod=zeros(2,2) #will end up being in transverse gauge
        @inbounds for i in 1:2
            @inbounds for j in 1:2
                prod[i,j] = sum(matrices[i,j]) #OUTER PRODUCT OF VECTORS    
            end
        end
        prod_TT=prod .- (prod[1,1]+prod[2,2]) .* [1 0; 0 1] #traceless
        return((factor/DL).*prod_TT)
    end
end;