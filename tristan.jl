### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 0cac30f1-9101-4ba3-accc-52621bc1d16f
begin
	#=import Pkg
    Pkg.add("PlutoTest")
	Pkg.add("PlutoUI")
	Pkg.add("PlutoTeachingTools")
	Pkg.add("PlutoLinks")=#
	using PlutoUI, PlutoTest, PlutoTeachingTools, PlutoLinks
	using ForwardDiff
	#package used for chain rules of large derivatives
	using OrdinaryDiffEq
	using SciMLBase
	#for integrating differential equations
	using PlutoTest
	include("victor-ftns_only.jl")
end;

# ╔═╡ ac043976-9ee2-4633-915b-1043a4771164


# ╔═╡ 647e477c-681f-43e7-9794-38d19bca6098


# ╔═╡ d06e9c81-4436-4054-a7cc-c85b7ee13df9
md"We model the neutron star as accreting conservatively until it reaches its Eddington limit, at which point it accretes at the Eddington limit (e.g., $\dot M _{\rm ns}=\min (-\dot M_{\rm WD}, \dot M_{\rm NS,\,Edd})$. However, as we shall see later, we require our functions to be differentiable. Thus, we use a 'smooth minimum' defined using the Boltzmann operator $S(x_1,x_2)=\frac{x_1e^{\alpha x_1}+x_2e^{\alpha x_2}}{e^{\alpha x_1}+e^{\alpha x_2}}$ for $\alpha<1$."

# ╔═╡ 07ca0c95-f9b3-42ad-8a49-3d305002b20f


# ╔═╡ 4f1118f3-e8f8-479f-819a-db201f36e447


# ╔═╡ 86648226-52f9-4951-bd85-bebc21018567
"""Smooth approximation to the minimum function"""
function Boltzmann_Max(x1::Float64, x2::Float64; α::Float64=-4.)::Float64
    return((x1*exp(α*x1)+x2*exp(α*x2))/(exp(α*x1)+exp(α*x2)))
end

# ╔═╡ 89014777-e6e0-406a-863c-72a310e37f7e
md"Testing for accuracy within 10% of true minimum value:"

# ╔═╡ 45d20cd8-a5e6-4709-b309-bfb3da22ed67
begin
	for i in 1:100
		test_x1, test_x2=rand(2)
		@test Boltzmann_Max(test_x1,test_x2) ≈ min(test_x1,test_x2) atol=0.1
	end
end

# ╔═╡ 4135feba-248d-4880-9e05-e3e85dae9963


# ╔═╡ 2fb18e7f-30eb-4dff-b30f-f5086ec3e02d
begin
	""" Smooth accretion rate of the neutron star """
	function dM_NS(a::Float64, R_WD::Float64, M_WD::Float64, M_NS::Float64;
			A::Float64=10, mu_e::Float64=2.0, eta_acc::Float64=0.1, α::Float64=4)::Float64
	    t_dM_WD=decretion_rate(a,a,R_WD,M_WD,M_Ns,A=A)
	    t_edd=eddington_rate(M_NS, mu_e=mu_e, eta_acc=eta_acc)
	    return(Boltzmann_Max(-t_dM_WD,t_edd,α=α))
	end

	""" Time derivative of angular parameter  θ"""
	function dθ(J::Float64,a::Float64,M_WD::Float64, M_NS::Float64)::Float64
		@assert M_WD > 0 # In production-ready code, consider commenting this out
		@assert M_NS > 0 # In production-ready code, consider commenting this out
		@assert r > 0 # In production-ready code, consider commenting this out
		return(J/(a*mu(M_WD,M_NS)))
	end   
end;

# ╔═╡ 4bbe25ad-bce2-4305-bd83-6dba2dc1f28c
md"Below, we define the quadrupole moment in terms of polar coordinates. (Note: As we are only concerned with higher derivatives of $I$ which are symmetric and carry 'plus/times' polarization, it is sufficient only to define $I_{xx/xy}$.) "

# ╔═╡ e90a2a13-e9af-4160-a307-cc7df36e280e
begin
	"""Given polar orbital parameters, calculates I_{xx} """
	function I_xx(M_WD::Float64, M_NS::Float64, r::Float64, θ::Float64)::Float64
	    @assert M_WD > 0 # In production-ready code, consider commenting this out
		@assert M_NS > 0 # In production-ready code, consider commenting this out
	    @assert r > 0 # In production-ready code, consider commenting this out
	    return (mu(M_WD, M_NS) * (x(r,θ))^2)
	end           

	"""Given polar orbital parameters, calculates I_{xy} """
	function I_xy(M_WD::Float64, M_NS::Float64, r::Float64, θ::Float64)::Float64
	    @assert M_WD > 0 # In production-ready code, consider commenting this out
		@assert M_NS > 0 # In production-ready code, consider commenting this out
	    @assert r > 0 # In production-ready code, consider commenting this out
	    return mu(M_WD, M_NS) * x(r,θ) * y(r,θ)
	end                
end;

# ╔═╡ 3fb94d9d-948f-4620-8621-e0f3d11e99cf


# ╔═╡ 84c2d206-599a-42e1-ba61-6bc06cb62f65


# ╔═╡ 7f29652e-29d7-43dc-8ad3-a1dae2f39991
md"To calculate the effect of the gravitational wave backreaction on the orbitm we must calculate derivatives up to $\dddot I^{ij}$. In principle, given the functions in the math section of README, it is possible to take these derivatives analytically by hand. In practice, this is a bit gruesome. We introduce autodifferentiation to take the necessary higher derivatives given that we can calculate the first time derivatives. For a function $f\equiv f(\vec{theta})$, we can use the chain rule to write: $\dot f =\sum_{i}\frac{\partial f}{\partial \theta^i}\dot\theta_i$
which is the basis for `time_derivative`. This gives that $\dot f\equiv \dot f(\vec{theta},\dot{\vec\theta})$, so to calculate $\ddot f$ involves taking the gradient over higher derivatives of the initial parameters. In general, $f^{(n)} =\sum_{j=0}^{n-1}\sum_{i}\frac{\partial f}{\partial (\theta^{i})^{(n)}}(\theta^{i})^{(n+1)}$ This is the basis for `next_derivative`.
"

# ╔═╡ 61544091-4d40-434b-bf52-de66d8c35c7a


# ╔═╡ 25998cc0-96e9-476c-8853-d102a58c0334


# ╔═╡ 0a57469d-418f-47b0-bf56-828d4bc529e1


# ╔═╡ 694fe531-e34e-4c27-94cb-ae3f1e14aa90
"""Given a function, its vectorized parameters, and the time derivatives of those parameters, returns the time derivative of the function as calculated via the chain rule."""
function time_derivative(ftn::Function, ind_vars::Vector{Float64}, derivatives::Vector{Float64})::Float64
    #ftn is the function to be taken a time derivative of
    #ind_vars is an array with the values of the arguments of ftn
    #derivatives is an array with the known time derivatives of ind_vars
    
    @assert length(ind_vars)==length(derivatives)
    
    grad=ForwardDiff.gradient(v->ftn(v[1],v[2],v[3],v[4]),ind_vars)
    #Packages ftn to take an array of parameters as variables (e.g., ind_vars)
    #then calculates the gradient of the function
    
    return(dot(grad,derivatives))#chain rule
end

# ╔═╡ 55048669-1a54-4ff9-acb3-194fd69dfc00


# ╔═╡ 621d6981-fed6-4665-85a2-f5df3a06519e
#=Higher derivatives=#
#=takes in values and time derivatives to calculate second derivative=#
"""Given a function, defines its next derivative based on the function's parameters and their time derivatives. `n` determines how many time derivatives of \vec{θ} (i.e., the parameters to be integrated) are parameters of the function itself. """
function next_derivative(in_ftn::Function; n::Integer=1)::Function
    #returns a function which is the next higher derivative of the in_ftn
    #n is the total derivatives past the given analytic form
    function ftn(gen_theta::Vector{Float64}, gen_dtheta::Vector{Float64})::Float64
        @assert length(gen_theta)==length(gen_dtheta)==6*n 
        @assert gen_theta[7:end]==gen_dtheta[1:end-6]
        #requires gen_theta to be [theta, dtheta, d...dtheta] for n=number of ds
        #requires gen_theta to be [dtheta, ddtheta d...ddtheta] for n+1=number of ds
        return(ftn(in_ftn, gentheta, gen_dtheta))
    end
    return ftn
end

# ╔═╡ 24cb207f-cd3a-4e50-ad02-6318c607738c
md"The following functions wrap the original derivatives in `victor.jl` in vectors to make life easier when taking derivatives due to mass transfer only. In particular, `theta` is of the form `[a, M_WD, M_NS, R_WD, J, θ]`."

# ╔═╡ 731205f1-64f9-4dd5-809f-a1cd20eabd47


# ╔═╡ 528d99bb-0800-41d7-812b-8685a18c2ec0


# ╔═╡ 3b211714-36b0-4b0d-9b62-6ae2939f66d9
begin
	#=The following functions wrap the original derivatives in vectors to make life easier when taking derivatives due to mass transfer only=#    
	function da_vector(theta::Vector{Float64}; 
			A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
			a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(a_dot(a,a,M_WD,M_NS,R_WD,
	            A=A,mu_e=mu_e,eta_acc=eta_acc,
	            a_ring_frac=a_ring_frac,mode=mode))
	end
	
	function dM_WD_vector(theta::Vector{Float64}; A::Float64=10.0)::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(decretion_rate(a, r, R_WD, M_WD, M_NS, A=A))
	end
	
	function dM_NS_vector(theta::Vector{Float64};
	        A::Float64=10, mu_e::Float64=2.0, eta_acc::Float64=0.1, α::Float64=4)::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(dM_NS(a, R_WD, M_WD, M_NS;
	        A=a, mu_e=mu_e, eta_acc=eta_acc, α=α))
	end
	
	function dR_WD_vector(theta::Vector{Float64}; mu_e::Float64=2.)::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return( R_WD(M_WD, mu_e=mu_e))
	end
	    
	function dJ_vector(theta::Vector{Float64}; 
	               A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
	               a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(J_dot(J, a, a, R_WD, M_WD, M_NS, 
	               A=A, mu_e=mu_e, eta_acc=eta_acc, 
	               a_ring_frac=a_ring_frac, mode=mode))
	end
	
	function dθ_vector(theta::Vector{Float64})::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(dθ(J,a,M_WD, M_NS))
	end
	
	
	#=Same treatment, but for quadrupole moments=#
	function Ixx_vector(theta::Vector{Float64})::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(I_xx(M_WD, M_NS, a, θ))
	end
	
	function Ixy_vector(theta::Vector{Float64})::Float64
	    a, M_WD, M_NS, R_WD, J, θ=theta
	    return(I_xy(M_WD, M_NS, a, θ))
	end
end

# ╔═╡ 8b8cf6b2-70ef-4c0c-b729-551637947f0c
md"The below block of code uses the vectorized derivatives above to output functions which take higher time derivatives of $\vec{\theta}$ and $I_{ij}$ as a function of $\theta$. "

# ╔═╡ 9878e729-0b7c-4684-a0fd-e0cf73465199
#=all necessary higher derivatives to calculate up to dddI=#
begin
    dda_vector=next_derivative(da_vector)
    ddM_WD_vector=next_derivative(dM_WD_vector)
    ddM_NS_vector=next_derivative(dM_WD_vector)
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

# ╔═╡ 0e80bd6c-80af-45b3-8637-22fa8bfe8a5c
#=bundling derivatives together which will prove convenient later=#
begin
    first_derivatives=[da_vector, dM_WD_vector, dM_WD_vector, dR_WD_vector, dJ_vector, dθ_vector]
    second_derivatives=[dda_vector, ddM_WD_vector, ddM_WD_vector, ddR_WD_vector, ddJ_vector, ddθ_vector]
    third_derivatives=[ddda_vector, dddM_WD_vector, dddM_WD_vector, dddR_WD_vector, dddJ_vector, dddθ_vector]
end;

# ╔═╡ 7f728824-950f-4631-bb6b-2b4654f6dc7f
#=Constructs a symmetric 2x2 tensor from a plus and times polarization=#
function tensor_plus_times(plus::Float64, times::Float64)::Vector{Float64}
    return([plus times;
            times -plus        
    ])
end

# ╔═╡ 5098b622-aacf-4df1-8341-637e8eac60e6
md"The below block of code is central to the interface between the mass transfer and radiation reaction effects. It also sits in the innermost loop of our integrator. This is a good place to peer review."

# ╔═╡ 6208496a-07a7-4f27-9994-98a3347f018b
#=Central Engine of the integrator=#
begin
	"""Calculates the second derivative of the quadrupole moment `t_ddI`  due to kinematics of keplerian motion and mass transfer `dtheta`. Also returns higher derivatives of theta useful for calculating the third derivative of the quadrupole moment."""
	function ddI_helper(theta::Vector{Float64})::Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Array{Float64,2}}
	    #calulates a list of values which will be helpful in the next two functions
	    dtheta=zeros(6)
	    #calculating first time derivatives due to MT only
	    for i in 1:6
	        dtheta[i]=first_derivatives(theta)
	    end
	
	    ddtheta=zeros(6)
	    #calculating second time derivatives due to MT only
	    for i in 1:6
	        ddtheta[i]=second_derivatives(theta,dtheta)
	    end
	
	    theta_dtheta=vcat(theta,dtheta)
	    dtheta_ddtheta=vcat(dtheta,ddtheta)
	
	    if not ddI_only
	        dddtheta=zeros(6)
	        #calculating third time derivatives due to MT only
	        for i in 1:6
	            dddtheta[i]=theta_derivatives(theta_dtheta,dtheta_ddtheta)
	        end
	    end
	
	
	    
	    #using orbital elements and first and second time derivatives 
	    #to calculate second derivative of quadrupole moment
	    #(for use in calculating RR acceleration)
	    
	    
	    t_ddI=tensor_plus_times(ddI_xx(theta_dtheta,dtheta_ddtheta),
	        ddI_xy(theta_dtheta,dtheta_ddtheta)
	    )
	    return(dtheta, ddtheta, theta_dtheta, dtheta_ddtheta, t_ddI)
	
	end

	"""Returns the second derivative of the quadrupole moment as calculated in `ddI_Helper`. """
	function ddI_from_theta(theta::Vector{Float64})::Array{Float64,2}
	    return(ddI_helper(theta)[5])
	end
	    
	
	"""Calculates the third derivative of the quadrupole moment and returns it along with the time derivative of theta and the second derivative of the quadrupole moment as calculated in `ddI_Helper`."""
	function ddI_dddI_from_theta(theta::Vector{Float64};ddI_only::Bool=false)::Tuple{Vector{Float64},Array{Float64,2},Array{Float64,2}}
	    #returns time derivatives of theta and second and third derivatives of 
	    #the quadrupole moment due to MT effects
	    
	    dtheta, ddtheta, theta_dtheta, dtheta_ddtheta, t_ddI = ddI_helper(theta)
	    #using orbital elements and first through third time derivatives  
	    #to calculate third derivative of quadrupole moment 
	    #(for use in calculating RR acceleration)
	    theta_ddtheta=vcat(theta_dtheta,ddtheta)
	    dtheta_dddtheta=vcat(dtheta_ddtheta,dddtheta)
	
	    t_dddI=tensor_plus_times(dddI_xx(theta_ddtheta,dddtheta),
	        dddI_xy(theta_ddtheta,dddtheta)
	    )
	    return(dtheta,t_ddI,t_dddI)
	end
end;

# ╔═╡ bc193a3e-052e-4860-854d-398fa08deeec


# ╔═╡ 120ff400-43a9-4be6-b1a8-2017ae0b111a
""" An "in place" implementation of the full time derivative of theta including effects from both mass transfer and radiation reaction.""" 
function acceleration!(du::Vector{Float64}, theta::Vector{Float64},p,t)
    dtheta,t_ddI,t_dddI=ddI_dddI_from_theta(theta)

    #forces due to RR
    t_dJ_rr=dJ_rr(t_ddI,t_dddI)
    t_da_rr=2.0 * theta[1] * dJ_rr(I2, I3) / theta[5]

    #adding RR forces to MT forces
    dtheta[1]+=t_da_rr
    dtheta[5]+=t_da_rr

    #mutate du
    du .= dtheta
    return(dtheta)
end;

# ╔═╡ 52b828d3-a6b7-4aa6-abe5-5c23ef848bb2


# ╔═╡ a63692e2-f652-4f35-b5c8-57ed131d676e
md"The next several blocks define functions intended to initiate the orbit beginning at mass transfer for initial masses, take an integration step forward calculating waveform, and integrate until the White Dwarf dissipates completely"

# ╔═╡ 72bfe0e8-11f2-47f6-97ba-09708af3fad3
#=The following functions are defined to initiate the orbit beginning at mass transfer for initial masses, take an integration step forward calculating waveform, and integrate until the White Dwarf dissipates completely=#

begin
	"""Calculates in km the orbital separation for the white dwarf to overflow its Roche Lobe and begin mass transfer"""
	function RL_contact(M_WD::Float64,M_NS::Float64;mu_e::Float64=2.)::Float64
	    #returns separation at which a donor WD overflows its
	    #Roche Lobe
	    R1=R_WD(M_WD,mu_e=mu_e)
	    RL_a1=Roche_Limit(1,M1,M2) #roche limit at separation of 1 km
	    return(R1 / RL_a1)#initial separation in km
	end

	"""Calculates the angular momentum of the system at initial Roche Lobe overflow """	
	function circular_J(a::Float64,M_WD::Float64,M_NS::Float64)::Float64
	    #returns angular momentum for a circular orbit
	    return((2 * a * G * M_WD ^ 2 * M_NS ^ 2 / (M1 + M2))^(1/2))
	end
	"""Defines theta at initial Roche Lobe overflow for  """
	function initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2.)::Float64
	    #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
	    #M_NS=1.35 M_\odot is typical NS mass
	    @assert M_WD<1.4 #Chandrasekar mass, beyond which we cannot calculate R_WD(see initial function)
	    a=RL_contact(M_WD,M_WD,mu_e=mu_e)
	    J=circular_J(a,M_WD,M_WD)
	    R1=R_WD(M_WD)
	    return([a, M_WD, M_NS, R_1, J, 0.])
	end
end

# ╔═╡ a7a7b50c-d555-4a1a-9197-963fac13aca5
md"Below is the heart of the integration process"

# ╔═╡ a8f954bf-b8d6-45f1-b836-d9a4c977f7db

begin
	"""Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken, the updated theta parameters, and the instantaneous second derivative of the quadupole moment at the end of the timestep."""
	function one_step(theta::Float64; dt_period::Float64=1e-2, dt::Float64=-1)::Tuple{Float64, Vector{Float64}, Array{Float64,2}}
	    #integrates forward one step with time interval either given explicitly (dt!=-1)
	    #or as a fraction of the orbital period determined by dt_period
	    if dt==-1
	        dt=dt_period*period(theta[2]+theta[3], theta[1])
	    end
	    prob = ODEProblem(acceleration,theta,(0,dt),[])
	    updated_theta = solve(prob, alg=KahanLi8()).u[end]
	
	    #=Calculating ddI at the end of the step=#
	        
	    t_ddI=ddI_from_theta(updated_theta)
	    return(dt,updated_theta,t_ddI)
	end

	"""Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. """
	function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
	        orbit_limit::Float64=1e5, step_limit::Int=10^7)::Tuple{Vector{Float64}, Array{Float64,3}}
	    #returns time series second derivatives of the quadrupole moment
	    t=0.
	    theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
	    ts=zeros(step_limit)
	    ddIs=zeros(step_limit,2,2)
	    i=0
	    time_limit=orbit_limit*period(M_WD+M_NS,theta[1])
	    
	    while t<time_limit && i<step_limit && theta[2]>0
	        #integrate until time or step limit are elapsed or WD completely dissipates
	        dt,t_theta,t_ddI=one_step(theta)
	        i+=1
	        t+=dt
	        ts[i]=t
	
	        theta=t_theta
	        ddIs[i, :, :].=ddI_from_theta(theta)
	    end
	
	    return(ts[1:i],ddIs[1:i, :, :])
	end
end;

# ╔═╡ 8f914e37-11f2-4620-b469-abc005e98312
integrate(1.0,step_limit=10)

# ╔═╡ 496e747c-7e5a-498b-ae70-604b9a64faef


# ╔═╡ 1f96a986-9dcb-4614-8720-44a193d5524e


# ╔═╡ bbee9e13-7a91-4662-ab1c-bd3ef245d1db


# ╔═╡ 43e3aa25-661c-40db-a638-f39c03f368b4


# ╔═╡ 7ec5ae9d-5a26-42f4-b1bd-17f61d9a3acf
md"Note: For peer review, we do NOT recommend running `integrate` fully (until total disruption) as this may take upward of 10 hours!"

# ╔═╡ af11b1f5-c481-4b76-a747-085d210af7d6
"""Given a luminosity distance `DL` in units of Mpc and an orientation parameterized by an inclination `i` and azimuthal angle `φ`, returns a function which calculates the instantaneous gravitational strain waveform based on the quadrupole formula."""
function gravitational_waveform_ftn(i::Float64, φ::Float64, DL::Float64 )::Function
    #fixes a function to return gravitational wave given orientation and distance
    #DL is in Mpc
    @assert 0<=i<=π
    @assert 0<=φ<2π
    @assert DL<10000 #Ensure correct units for DL
    e1=[cos(i)*cos(φ), sin(i)*sin(φ), cos(i)]
    e2=[-sin(φ),  cos(φ),  0]
    es=[e1,e2]

    matrices=zeros(4,2,2)
    for i in 1:2
        for j in 1:2
            matrices[2*(i-1)+j ,:,:] = es[i] * es[j] #OUTER PRODUCT OF VECTORS
        end
    end

    factor= 2*G*3.2e-20/c^4 #for conversion from Mpc to km
    
    function ftn(ddI::Vector{Float64})::Vector{Float64}
        #given fixed outer products, calculate strain
        prod=zeros(2,2) #will end up being in transverse gauge
        for i in 1:2
            for j in 1:2
                prod[i,j] = sum(matrices[i,j]) #OUTER PRODUCT OF VECTORS    
            end
        end
        prod_TT=prod .- (prod[1,1]+prod[2,2]) .* [1 0; 0 1] #traceless
        return((factor/DL).*prod_TT)
    end
end;

# ╔═╡ 37b97235-ee25-45fd-a766-312debd46bd8


# ╔═╡ c8bae4b5-f4e6-429c-b1cd-4a2120eb3680


# ╔═╡ 9a3e04f2-bc3e-4f1e-80ec-4106635ecb01


# ╔═╡ 0797c258-aeec-4154-aee8-5b0359e6d0d1
#=
function dI_xx(M_WD::Float64, M_NS::Float64, a::Float64, θ::Float64,
        t_dM_WD::Float64,t_dM_NS::Float64,t_da::Float64,t_dθ::Float64)::Float64
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert a > 0 # In production-ready code, consider commenting this out
            
    grad=ForwardDiff.gradient(v->I_xx(v[1],v[2],v[3],v[4]), 
        [M_WD, M_NS, a, θ])#gradient of I_xx  wrt tracked parameters

    return(dot(grad,[ t_dM_WD,t_dM_NS,t_da,t_dθ])) 
end

function dI_xy(M_WD::Float64, M_NS::Float64, a::Float64, θ::Float64,
        t_dM_WD::Float64,t_dM::Float64_NS,t_da::Float64,t_dθ::Float64)::Float64
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert a > 0 # In production-ready code, consider commenting this out
            
    grad=ForwardDiff.gradient(v->I_xy(v[1],v[2],v[3],v[4]), 
        [M_WD, M_NS, a, θ])#gradient of I_xx  wrt tracked parameters

    return(dot(grad,[ t_dM_WD,t_dM_NS,t_da,t_dθ])) 
end
=#

                
#=function dI_xx(M_WD::Float64, M_NS::Float64, a::Float64, θ::Float64, R_WD::Float64, J::Float64; 
                    A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
               a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert a > 0 # In production-ready code, consider commenting this out
            
    grad=ForwardDiff.gradient(v->I_xx(v[1],v[2],v[3],v[4]), 
        [M_WD, M_NS, a, θ])#gradient of I_xx  wrt tracked parameters

    #time derivatives of  M_WD,M_NS,a,θ       
    t_dM_WD=decretion_rate(a, a, R_WD, M_WD, M_NS, A=A)
    t_dM_NS=dM_NS(a,R_WD,M_WD,M_NS)
    t_da=a_dot(a,a, R_WD, M_WD, M_NS, 
               A=A, mu_e:=mu_e, eta_acc=eta_acc, a_ring_frac=1.0, mode="Isotropic")
    t_dθ=dθ(J,a,M_WD,M_NS)

    return(dot(grad,[ t_dM_WD,t_dM_NS,t_da,t_dθ]))                    
            
function dI_yy(M_WD::Float64, M_NS::Float64, a::Float64, θ::Float64, R_WD::Float64, J::Float64; 
                    A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
               a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert a > 0 # In production-ready code, consider commenting this out
            
    grad=ForwardDiff.gradient(v->I_yy(v[1],v[2],v[3],v[4]), 
        [M_WD, M_NS, a, θ])#gradient of I_yy  wrt tracked parameters

    #time derivatives of  M_WD,M_NS,a,θ       
    t_dM_WD=decretion_rate(a, a, R_WD, M_WD, M_NS, A=A)
    t_dM_NS=dM_NS(a,R_WD,M_WD,M_NS)
    t_da=a_dot(a,a, R_WD, M_WD, M_NS, 
               A=A, mu_e:=mu_e, eta_acc=eta_acc, a_ring_frac=1.0, mode="Isotropic")
    t_dθ=dθ(J,a,M_WD,M_NS)

    return(dot(grad,[ t_dM_WD,t_dM_NS,t_da,t_dθ])) =# 

# ╔═╡ Cell order:
# ╠═0cac30f1-9101-4ba3-accc-52621bc1d16f
# ╠═ac043976-9ee2-4633-915b-1043a4771164
# ╠═647e477c-681f-43e7-9794-38d19bca6098
# ╠═d06e9c81-4436-4054-a7cc-c85b7ee13df9
# ╠═07ca0c95-f9b3-42ad-8a49-3d305002b20f
# ╠═4f1118f3-e8f8-479f-819a-db201f36e447
# ╠═86648226-52f9-4951-bd85-bebc21018567
# ╠═89014777-e6e0-406a-863c-72a310e37f7e
# ╠═45d20cd8-a5e6-4709-b309-bfb3da22ed67
# ╠═4135feba-248d-4880-9e05-e3e85dae9963
# ╠═2fb18e7f-30eb-4dff-b30f-f5086ec3e02d
# ╠═4bbe25ad-bce2-4305-bd83-6dba2dc1f28c
# ╠═e90a2a13-e9af-4160-a307-cc7df36e280e
# ╠═3fb94d9d-948f-4620-8621-e0f3d11e99cf
# ╠═84c2d206-599a-42e1-ba61-6bc06cb62f65
# ╟─7f29652e-29d7-43dc-8ad3-a1dae2f39991
# ╠═61544091-4d40-434b-bf52-de66d8c35c7a
# ╠═25998cc0-96e9-476c-8853-d102a58c0334
# ╠═0a57469d-418f-47b0-bf56-828d4bc529e1
# ╠═694fe531-e34e-4c27-94cb-ae3f1e14aa90
# ╠═55048669-1a54-4ff9-acb3-194fd69dfc00
# ╠═621d6981-fed6-4665-85a2-f5df3a06519e
# ╠═24cb207f-cd3a-4e50-ad02-6318c607738c
# ╠═731205f1-64f9-4dd5-809f-a1cd20eabd47
# ╠═528d99bb-0800-41d7-812b-8685a18c2ec0
# ╠═3b211714-36b0-4b0d-9b62-6ae2939f66d9
# ╠═8b8cf6b2-70ef-4c0c-b729-551637947f0c
# ╠═9878e729-0b7c-4684-a0fd-e0cf73465199
# ╠═0e80bd6c-80af-45b3-8637-22fa8bfe8a5c
# ╠═7f728824-950f-4631-bb6b-2b4654f6dc7f
# ╠═5098b622-aacf-4df1-8341-637e8eac60e6
# ╠═6208496a-07a7-4f27-9994-98a3347f018b
# ╠═bc193a3e-052e-4860-854d-398fa08deeec
# ╠═120ff400-43a9-4be6-b1a8-2017ae0b111a
# ╠═52b828d3-a6b7-4aa6-abe5-5c23ef848bb2
# ╠═a63692e2-f652-4f35-b5c8-57ed131d676e
# ╠═72bfe0e8-11f2-47f6-97ba-09708af3fad3
# ╠═a7a7b50c-d555-4a1a-9197-963fac13aca5
# ╠═a8f954bf-b8d6-45f1-b836-d9a4c977f7db
# ╠═8f914e37-11f2-4620-b469-abc005e98312
# ╠═496e747c-7e5a-498b-ae70-604b9a64faef
# ╠═1f96a986-9dcb-4614-8720-44a193d5524e
# ╠═bbee9e13-7a91-4662-ab1c-bd3ef245d1db
# ╠═43e3aa25-661c-40db-a638-f39c03f368b4
# ╠═7ec5ae9d-5a26-42f4-b1bd-17f61d9a3acf
# ╠═af11b1f5-c481-4b76-a747-085d210af7d6
# ╠═37b97235-ee25-45fd-a766-312debd46bd8
# ╠═c8bae4b5-f4e6-429c-b1cd-4a2120eb3680
# ╠═9a3e04f2-bc3e-4f1e-80ec-4106635ecb01
# ╟─0797c258-aeec-4154-aee8-5b0359e6d0d1
