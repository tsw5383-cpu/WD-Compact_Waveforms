### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 0cac30f1-9101-4ba3-accc-52621bc1d16f
#=begin
	#=import Pkg
    Pkg.add("PlutoTest")
	Pkg.add("PlutoUI")
	Pkg.add("PlutoTeachingTools")
	Pkg.add("PlutoLinks")
	Pkg.add("OrdinaryDiffEq")
	Pkg.add("SciMLBase")
	Pkg.add("ForwardDiff")
	Pkg.add("BenchmarkTools")=#

	#=import Pkg
	=#

	#Things I need for test macro
	using LinearAlgebra
	using BenchmarkTools
	using Profile, ProfileCanvas 
	using PlutoUI, PlutoTeachingTools, PlutoLinks

	
	using ForwardDiff
	#package used for chain rules of large derivatives
	
	using OrdinaryDiffEq
	using SciMLBase
	#for integrating differential equations

	
	
end;=#

begin
	using LinearAlgebra
	using BenchmarkTools
	using Profile, ProfileCanvas 
	using PlutoUI, PlutoTeachingTools, PlutoLinks
	using Test

	using ForwardDiff, OrdinaryDiffEq, SciMLBase
	#packages used for chain rules of large derivatives
	#and integrating differential equations

	using CSV, DataFrames
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

# ╔═╡ 6c59b881-5668-4894-950b-2c7bee5367ef
md"Importing packages is causing me a world of problems."

# ╔═╡ 511606a0-8594-45a3-af99-29a3ffbebe68


# ╔═╡ 8edaf48b-d27d-43c9-8a0c-e149654b6cd0


# ╔═╡ c79771dc-4bf2-4c03-9f84-dd6532c98bfd


# ╔═╡ 33f9cc9a-a207-4a09-bfa3-0f54546c05c3


# ╔═╡ e2833dab-01a3-4d20-a42d-7d182f9e610d


# ╔═╡ 1db7cf48-6630-4701-b163-c2fcc632a050


# ╔═╡ 5519a048-4a4a-4091-af9e-8b2dc444775b


# ╔═╡ c13a77be-813e-413f-b0b2-23cf7b5f309a


# ╔═╡ 14fedf73-b9ca-4a37-a86c-e6f4c0a59ece


# ╔═╡ abca47ab-091a-4677-8d5c-b13572bb7baa


# ╔═╡ b97703e5-e4f3-46a9-8819-81568e8582d6


# ╔═╡ 4bef6030-10c5-4497-9ec2-edabee8432d1


# ╔═╡ 17be5889-9440-4adb-965f-17de4875a922


# ╔═╡ ce80bb70-d60e-4bdf-8f7c-16421fd90138


# ╔═╡ fa5a36f9-b1e2-4461-bb18-a3172f4b7af7


# ╔═╡ 421df964-fd81-471d-bff5-570f3f409f53


# ╔═╡ 2ea8a9da-f653-41f8-97df-f2909f7ce1a2


# ╔═╡ 25d4f52f-7f37-41c9-9230-37dfb0e71f22


# ╔═╡ 104e6849-3bae-467a-8630-57297b116314


# ╔═╡ 3ef782fd-4d6a-4c66-a976-5e995e2d703d
@__DIR__

# ╔═╡ 13bef496-5bb6-4dbe-a6e4-d18e664dfc55


# ╔═╡ 9a2f1c1e-46ca-41b8-9d11-d165bb86a918


# ╔═╡ 3e28aef1-3662-47b6-8d59-9a1c97313d5b


# ╔═╡ 190a4e03-541a-4857-879d-7093afdfddc2


# ╔═╡ d3985c1e-dde5-4700-b948-b6fd04c80819


# ╔═╡ 40d859c3-a08e-49a2-8b8b-34ad7a0a6fa3


# ╔═╡ 35d54cfa-7d0a-4a48-94bc-a45327b0abcc


# ╔═╡ 1620f45c-a715-432c-8a2d-f673bb022e2c


# ╔═╡ 96c169cc-f4d7-4885-ab1f-a692f5f39f91


# ╔═╡ d0cc1fbf-6828-4bf5-8bb6-1f9d97161a24


# ╔═╡ 1ffcf92d-958a-44a2-88cc-1c7ab262925b


# ╔═╡ 8aca584e-0052-4e4d-9e30-823204f9595e


# ╔═╡ fc8da2c6-e7b8-47cd-af05-b13e74909dda


# ╔═╡ e4452b17-0722-4ea0-804b-7f2797e17a82


# ╔═╡ ac043976-9ee2-4633-915b-1043a4771164


# ╔═╡ 647e477c-681f-43e7-9794-38d19bca6098


# ╔═╡ d06e9c81-4436-4054-a7cc-c85b7ee13df9
md"We model the neutron star as accreting conservatively until it reaches its Eddington limit, at which point it accretes at the Eddington limit (e.g., $\dot M _{\rm ns}=\min (-\dot M_{\rm WD}, \dot M_{\rm NS,\,Edd})$. However, as we shall see later, we require our functions to be differentiable. Thus, we use a 'smooth minimum' defined using the Boltzmann operator $S(x_1,x_2)=\frac{x_1e^{\alpha x_1}+x_2e^{\alpha x_2}}{e^{\alpha x_1}+e^{\alpha x_2}}$ for $\alpha<0$."

# ╔═╡ 07ca0c95-f9b3-42ad-8a49-3d305002b20f


# ╔═╡ 4f1118f3-e8f8-479f-819a-db201f36e447


# ╔═╡ 86420e50-5f10-43aa-ac65-a8708be68dc6
"""Smooth approximation to the minimum function"""
function Boltzmann_Max(x1::T, x2::T; α=-10.) where T<:Real
    return((x1*exp(α*x1)+x2*exp(α*x2))/(exp(α*x1)+exp(α*x2)))
end

# ╔═╡ 6d4270fb-21e4-4e8c-b7e4-0e7bc83870c6


# ╔═╡ 94098f26-fc8d-4df5-84fe-9ebb7b8c994a
#=function Boltzmann_Min(x1, x2; α=-10.)
    Δ = α * (x2 - x1)
	out=0
    if Δ > 50
        out=x1           # exp(Δ) overflows → x1 << x2
    elseif Δ < -50
        out=x2           # exp(-Δ) overflows → x2 << x1
    else
        w = 1 / (1 + exp(Δ))
        out= w * x1 + (1 - w) * x2
	end		
	return out 
end=#

# ╔═╡ 256492bf-b826-49aa-b051-55991228703c


# ╔═╡ a3e1afb5-0a7b-41aa-b50f-dc959a36a1e6


# ╔═╡ 69081296-10a7-4628-a62d-b116548bdbae


# ╔═╡ bfbaf692-db1f-41f9-a147-95ccb8de9e1e


# ╔═╡ 1de07299-8623-49fe-bba3-136dc45e7c82


# ╔═╡ ea61a15f-1982-4469-8e0b-22f309d1bef5


# ╔═╡ 89014777-e6e0-406a-863c-72a310e37f7e
md"Testing for accuracy within 10% of true minimum value:"

# ╔═╡ 45d20cd8-a5e6-4709-b309-bfb3da22ed67
begin
	for i in 1:100
		test_x1, test_x2=rand(2)
		@test Boltzmann_Max(test_x1,test_x2) ≈ min(test_x1,test_x2) atol=0.1
	end
end

# ╔═╡ b2dcf719-078c-4ec0-bde6-243972cfc957


# ╔═╡ 4135feba-248d-4880-9e05-e3e85dae9963


# ╔═╡ 2fb18e7f-30eb-4dff-b30f-f5086ec3e02d
begin
	""" Smooth accretion rate of the neutron star """
	#=function dM_NS(a::T, R_WD::T, M_WD::T, M_NS::T;
			A=10., mu_e=2.0, eta_acc=0.1, α=-10.) where T<:Real
	    t_dM_WD=decretion_rate(a,a,R_WD,M_WD,M_NS,A=A)
		
	    t_edd=eddington_rate(M_NS, mu_e=mu_e, eta_acc=eta_acc)

		out=Boltzmann_Max(-1*t_dM_WD,t_edd,α=α)

		println(t_dM_WD)
		println(t_edd)
		println(out)
		@assert out>0
		@assert out<-1.5*t_dM_WD
		@assert out<1.5*t_edd
	    return(out)
	end=#

	function dM_NS(a::T, R_WD::T, M_WD::T, M_NS::T;
			A=10., mu_e=2.0, eta_acc=0.1, α=-10.) where T<:Real
	    t_dM_WD=decretion_rate(a,a,R_WD,M_WD,M_NS,A=A)
		
	    t_edd=eddington_rate(M_NS, mu_e=mu_e, eta_acc=eta_acc)
		return(min(-t_dM_WD,t_edd))
	end
	

	println(dM_NS(14088.7,5338.49,1.,1.35))

	"""Reduced mass in effective one body formalism"""
	function mu(M_WD::T, M_NS::T) where T<:Real
		return((M_WD+M_NS)/(M_WD+M_NS))
	end

	
	""" Time derivative of angular parameter  θ"""
	function dθ(J::T,a::T,M_WD::T, M_NS::T) where T<:Real
		@assert M_WD > 0 # In production-ready code, consider commenting this out
		@assert M_NS > 0 # In production-ready code, consider commenting this out
		@assert a > 0 # In production-ready code, consider commenting this out
		return(J/(a^2*mu(M_WD,M_NS)))
	end   
end;

# ╔═╡ f784cf7a-6cf3-48b0-b705-acf708b75db8


# ╔═╡ d43a93a0-3f71-4d7f-b548-6057a1513f87


# ╔═╡ d45576f6-f2ad-47e5-a4a7-6bf6f95d42c9


# ╔═╡ 557daf12-69f8-4f4f-894f-6bc435fa91ab


# ╔═╡ df23f515-ffc2-4965-b695-c3875c5bf29d


# ╔═╡ 760d0fda-81d0-4d36-9e3d-46912d0ddfa4


# ╔═╡ 6ac10460-5123-487d-a55d-6183de743a5e


# ╔═╡ d93a9205-0f8c-4b47-af4e-5ad0b6437353


# ╔═╡ 344fdbd0-96e3-4a34-8d9f-3b6c34cee2b3
dM_NS(1,1,1,1)

# ╔═╡ 5512bb36-6de6-41e6-93a4-d736f8b69b3d
dM_NS(14088.7,5338.49,1.,1.35)

# ╔═╡ a1e4b8ec-de2e-42ab-8d52-2f62cf5d92dd
dot_R_WD()

# ╔═╡ c16514a9-c297-451b-9e6d-1d21d15d0f28


# ╔═╡ d9e0ca74-8f4b-4d10-8255-7a1dae0d49d7


# ╔═╡ 7fc2e8cf-556f-421f-a993-34e845c4ed3f


# ╔═╡ 0757c720-c651-44bc-aca6-20798890e01b


# ╔═╡ 5a5761f7-eb8e-419f-8bb6-15780de928d2


# ╔═╡ 8f6900c9-352f-4deb-82e2-1bd99faffda4


# ╔═╡ bb73d75d-eeeb-4b38-8db2-5f9cd4a440c3


# ╔═╡ 15095435-54c0-4d8d-91c0-c278edd57700


# ╔═╡ 870a48df-44ba-4585-b422-c92a84737a12


# ╔═╡ 2e25bf2a-d4c0-41c6-b4f3-1fff15bd6032


# ╔═╡ 7ffbc67d-2ad6-4ca2-8e52-7bcd81089da1


# ╔═╡ a18ec6f3-6036-4e0e-b8f2-6156e0497e8c


# ╔═╡ 69640bb5-07a4-4257-9eea-6f16dfcce74c
dM_NS(1.01, 1.01, 1.01, 1.01)

# ╔═╡ 7faf665d-0bdc-499d-8cc6-72ebd7da30b7


# ╔═╡ c42a4b21-75f4-4246-a34c-c51fdd752d0a


# ╔═╡ 94e9ebac-2ffa-4e76-be6c-273f402064ae


# ╔═╡ 57900831-e4f4-4557-bacc-422af3cc2d2b


# ╔═╡ 02584e20-271e-4263-a0b9-33889bc18ae6


# ╔═╡ dd2c1373-c195-49b1-8f00-18894900bd85


# ╔═╡ c69bd40e-49a0-4f5c-b516-6d5f774dd75d


# ╔═╡ 2efd9365-9bcf-41ea-9e77-9346bf498bde


# ╔═╡ 6b59b741-3cfa-4db2-afce-e6f23ab1a06a


# ╔═╡ 4bbe25ad-bce2-4305-bd83-6dba2dc1f28c
md"Below, we define the quadrupole moment in terms of polar coordinates. (Note: As we are only concerned with higher derivatives of $I$ which are symmetric and carry 'plus/times' polarization, it is sufficient only to define $I_{xx/xy}$.) "

# ╔═╡ e90a2a13-e9af-4160-a307-cc7df36e280e
begin
	"""Given polar orbital parameters, calculates I_{xx} """
	function I_xx(M_WD::T, M_NS::T, r::T, θ::T) where T<:Real
	    @assert M_WD > 0 # In production-ready code, consider commenting this out
		@assert M_NS > 0 # In production-ready code, consider commenting this out
	    @assert r > 0 # In production-ready code, consider commenting this out
	    return (mu(M_WD, M_NS) * (x(r,θ))^2)
	end           

	"""Given polar orbital parameters, calculates I_{xy} """
	function I_xy(M_WD::T, M_NS::T, r::T, θ::T) where T<:Real
	    @assert M_WD > 0 # In production-ready code, consider commenting this out
		@assert M_NS > 0 # In production-ready code, consider commenting this out
	    @assert r > 0 # In production-ready code, consider commenting this out
	    return mu(M_WD, M_NS) * x(r,θ) * y(r,θ)
	end                
end;

# ╔═╡ aad547b0-ff39-4c6c-bd40-d47d132e3f4d


# ╔═╡ d877e869-f0ec-4b02-9bdc-8f8d6d4e1c6c


# ╔═╡ deb8fbb9-e4f9-44e3-87cf-e226319b8a2e


# ╔═╡ 3fb94d9d-948f-4620-8621-e0f3d11e99cf


# ╔═╡ 84c2d206-599a-42e1-ba61-6bc06cb62f65


# ╔═╡ 7f29652e-29d7-43dc-8ad3-a1dae2f39991
md"To calculate the effect of the gravitational wave backreaction on the orbitm we must calculate derivatives up to $\dddot I^{ij}$. In principle, given the functions in the math section of README, it is possible to take these derivatives analytically by hand. In practice, this is a bit gruesome. We introduce autodifferentiation to take the necessary higher derivatives given that we can calculate the first time derivatives. For a function $f\equiv f(\vec{\theta})$, we can use the chain rule to write: $\dot f =\sum_{i}\frac{\partial f}{\partial \theta^i}\dot\theta_i$
which is the basis for `time_derivative`. This gives that $\dot f\equiv \dot f(\vec{theta},\dot{\vec\theta})$, so to calculate $\ddot f$ involves taking the gradient over higher derivatives of the initial parameters. In general, $f^{(n)} =\sum_{j=0}^{n-1}\sum_{i}\frac{\partial f}{\partial (\theta^{i})^{(n)}}(\theta^{i})^{(n+1)}$ This is the basis for `next_derivative`.
"

# ╔═╡ 69e12b46-10d9-497b-9d2c-7b54e1658d7f


# ╔═╡ 61544091-4d40-434b-bf52-de66d8c35c7a


# ╔═╡ 25998cc0-96e9-476c-8853-d102a58c0334


# ╔═╡ 0a57469d-418f-47b0-bf56-828d4bc529e1
"""Fixes the higher derivatives of the function so that it is formatted in a way that ForwardDiff is comfortable with."""
function fix_higher_derivatives(ftn::Function,higher_derivatives::AbstractVector{T}) where {T} 
	function inner_ftn(gen_theta::AbstractVector{T}) where {T}
		return(ftn(gen_theta,higher_derivatives))
	end
	return ftn 
end;

# ╔═╡ 694fe531-e34e-4c27-94cb-ae3f1e14aa90
"""Given a function, its vectorized parameters, and the time derivatives of those parameters, returns the time derivative of the function as calculated via the chain rule."""
#=
function time_derivative_two_args(ftn::Function, ind_vars::AbstractVector{T}, derivatives::AbstractVector{T}) where {T}
    #ftn is the function to be taken a time derivative of
    #ind_vars is an array with the values of the arguments of ftn
    #derivatives is an array with the known time derivatives of ind_vars
	
    @assert length(ind_vars)==length(derivatives)
    
    grad=ForwardDiff.gradient(ftn,ind_vars)
    #Packages ftn to take an array of parameters as variables (e.g., ind_vars)
    #then calculates the gradient of the function
    
    return(dot(grad,derivatives))#chain rule
end=#

# ╔═╡ 80ee41e0-5978-47fe-8554-16477fb317d6


# ╔═╡ ba3f4fbe-1eb3-4040-8616-28f926bb6243


# ╔═╡ fdd55e22-89d2-465a-80b8-5c59f5f7129b


# ╔═╡ 28bdbf35-455c-4e01-a1bb-05c862abe361
"""Wraps the vectorized parameters and their derivatives into one input vector ind_vars_and_derivatives"""
function time_derivative(ftn::Function, ind_vars_and_derivatives::AbstractVector{T}) where {T}
    #ftn is the function to be taken a time derivative of
    #ind_vars is an array with the values of the arguments of ftn
    #derivatives is an array with the known time derivatives of ind_vars
	#print(length(ind_vars_and_derivatives))
    
    @assert mod(length(ind_vars_and_derivatives),6)==0

	ind_vars=ind_vars_and_derivatives[1:end-6]
	derivatives=ind_vars_and_derivatives[7:length(ind_vars_and_derivatives)]

	grad=ForwardDiff.gradient(ftn,ind_vars)
    
    return(sum(grad .* derivatives))#chain rule
end

# ╔═╡ 72720c31-b465-4755-b070-9f34311a0d97


# ╔═╡ e1591006-147d-484b-8559-763e239af504


# ╔═╡ 57002625-8ba1-43a1-b2e3-0094fa53c841


# ╔═╡ 5e8ec2f8-d73a-4a81-86dc-bb9c3b538caa


# ╔═╡ c1636ad3-16b4-43d1-b2bb-993885d11818


# ╔═╡ 44ba797c-48c8-45e2-87af-a6bf102066a3


# ╔═╡ 951f256e-c418-4f39-b4aa-0f4080caae8b


# ╔═╡ a238a819-a6ef-479c-912a-2ddc5812ce66


# ╔═╡ 72dd65d6-0b35-4d2e-a359-3abe9aa39c5f
begin
	"""A very simple (effectively two argument) function serving as a unit test for `time_derivative`
	"""
	function test_ftn(xs::AbstractVector{T}) where {T}
		x,y, x1, x2, x3, x4 =xs
		return(x+y^2)
	end;
	println(length([1., 2., 0. , 0., 0., 0., 1., 1., 0. , 0., 0., 0.]))
	println(time_derivative(test_ftn, [1., 2., 0. , 0., 0., 0., 1., 1., 0. , 0., 0., 0.]) )
	@test time_derivative(test_ftn, [1., 2., 0. , 0., 0., 0., 1., 1., 0. , 0., 0., 0.]) ≈ 5 atol=1e-15
end;

# ╔═╡ 92a94cfd-a3ff-4d62-af3b-dc8b0b57b0ba
methods(test_ftn)

# ╔═╡ 4606dfbb-e601-4da5-8ed0-e7cf1e5be3db


# ╔═╡ bd182ff6-48cc-445e-b067-e387b88ff08b


# ╔═╡ 4488d59c-5ef0-4d3e-8695-161bfefaf51c


# ╔═╡ 6951b817-887e-46b7-9e5c-e4840ffd6cf4


# ╔═╡ d71448b4-705e-468a-91b5-d3c87b599316


# ╔═╡ a3315e6b-846c-4f8b-96bf-aa2519d44ae7


# ╔═╡ 4a299a1e-afeb-4c7c-82a9-a3140a093719


# ╔═╡ 1efcd158-2476-43ca-9474-66ec7e391aa0


# ╔═╡ f52559d9-35d7-4464-bd88-26419fb5354f


# ╔═╡ d95b1547-717e-4ba4-8a1c-edea66926550


# ╔═╡ 11d2af20-64fd-41de-b33c-76122c047bcb


# ╔═╡ 621d6981-fed6-4665-85a2-f5df3a06519e
#=Higher derivatives=#
#=takes in values and time derivatives to calculate second derivative=#
"""Given a function, defines its next derivative based on the function's parameters and their time derivatives. `n` determines how many time derivatives of \vec{θ} (i.e., the parameters to be integrated) are parameters of the function itself. """
#=function next_derivative(in_ftn::Function; n::Integer=1)::Function
    #returns a function which is the next higher derivative of the in_ftn
    #n is the total derivatives past the given analytic form
    function ftn(gen_theta::AbstractVector{T}, gen_dtheta::AbstractVector{T})::Float64 where {T}
		print(in_ftn)
		println(length(gen_theta),length(gen_dtheta),6*n)
        @assert length(gen_theta)==length(gen_dtheta)==6*n 
        @assert gen_theta[7:end]==gen_dtheta[1:end-6]
		ind_vars_and_derivatives
        #requires gen_theta to be [theta, dtheta, d...dtheta] for n=number of ds
        #requires gen_theta to be [dtheta, ddtheta d...ddtheta] for n+1=number of ds
        return(time_derivative(in_ftn, gen_theta, gen_dtheta))
    end
    return ftn
end=#

# ╔═╡ 05576550-3ffb-4d2e-b89e-e4d00ed54d71
"""Given a function, defines its next derivative based on the function's parameters and their time derivatives. `n` determines how many time derivatives of \vec{θ} (i.e., the parameters to be integrated) are parameters of the function itself. """
function next_derivative(in_ftn::Function; n::Integer=1)
    #returns a function which is the next higher derivative of the in_ftn
    #n is the total derivatives past the given analytic form
	@assert n>0
    function ftn(gen_theta_and_derivatives::AbstractVector{T}) where {T}
		#print(in_ftn)
        @assert length(gen_theta_and_derivatives)==6*(n+1) 

		#outdated code
		#=gen_theta=gen_theta_and_derivatives[1:end-6]
		gen_dtheta=gen_theta_and_derivatives[end-5:end]=#
        #requires gen_theta to be [theta, dtheta, d...dtheta] for n=number of ds
        #requires gen_theta to be [dtheta, ddtheta d...ddtheta] for n+1=number of ds
        return(time_derivative(in_ftn, gen_theta_and_derivatives))
    end
    return ftn
end

# ╔═╡ 509ec580-b9b1-42dd-ab81-b72d1e1f7293


# ╔═╡ 19be5fd8-ea30-4cb7-93a4-cc2171651194


# ╔═╡ 5bdbc460-4df3-4d64-87d2-7820997d2317


# ╔═╡ b889ebe5-e259-4208-8685-cc32df37f58c


# ╔═╡ 95d0ce5f-fe9a-4cd2-bb0f-90eb4a5e1f39


# ╔═╡ d41d2eb7-701e-4893-8116-7d477effeef4


# ╔═╡ de0d4bc7-e2b8-4812-8b28-43bda221b17f


# ╔═╡ 31736097-efce-4c31-85ab-d657cb97822a


# ╔═╡ 19a2498d-e8c1-4394-bc7c-9dd2e6033b16


# ╔═╡ 96a595df-5e71-4d63-800f-f600ef4b591a


# ╔═╡ e681a7b9-9996-4105-ac4c-f1bd880f8235


# ╔═╡ 7d4b6bb6-2fcc-4669-9c4d-369b10b2644b


# ╔═╡ 7fcacf8f-b82f-48f0-87ab-9f11b7ac936f


# ╔═╡ f9e6f40b-ebb8-4711-99b2-6852b6a9c223


# ╔═╡ ec453c88-b887-472e-a98c-f96c3b896f08


# ╔═╡ 9d610030-0721-440d-a7dd-ea6007249228


# ╔═╡ 24cb207f-cd3a-4e50-ad02-6318c607738c
md"The following functions wrap the original derivatives in `victor.jl` in vectors to make life easier when taking derivatives due to mass transfer only. In particular, `theta` is of the form `[a, M_WD, M_NS, R_WD, J, θ]`."

# ╔═╡ 731205f1-64f9-4dd5-809f-a1cd20eabd47


# ╔═╡ 528d99bb-0800-41d7-812b-8685a18c2ec0


# ╔═╡ 3b211714-36b0-4b0d-9b62-6ae2939f66d9
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

		
	    return(a_dot(a,a,M_WD,M_NS,R_WD;
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

# ╔═╡ bf13e120-01a9-4026-90ee-29576bfa784c


# ╔═╡ 680002e1-a8e8-4fec-a37c-e6c92d8a6d2e


# ╔═╡ 59b217a8-0b07-442e-bcd4-7fda42d5d03a


# ╔═╡ d0db0f15-0632-4f3b-bdc0-af1342ae03a7


# ╔═╡ 2f9e1da8-3522-4aad-af07-00d13602b610


# ╔═╡ fd1ffa20-5780-4ee5-8f6d-248c02a1656e


# ╔═╡ 736914e5-11e5-40bb-8792-236e132699d5


# ╔═╡ 47c8e194-843f-4648-95c7-ad08a96cba19


# ╔═╡ 42d106c9-1b6f-42ed-966b-4c5fb40b7f5f


# ╔═╡ 2f6194b0-615f-4342-a49a-7364756b4d01


# ╔═╡ cf21fed8-294e-40ee-b29a-b1e0b673b834


# ╔═╡ a00e5b5a-488e-413e-aea1-f404dec591b5


# ╔═╡ 4ba2c0d2-c006-4a2d-ae05-48fc40948cf6


# ╔═╡ e2064499-24de-4c0a-a975-289dcf61b361


# ╔═╡ ecce3a10-b519-49d9-ab59-bd33926b5e3d


# ╔═╡ ad39dbb0-a043-4ed9-ada1-5a59612f57ec


# ╔═╡ e70e5b74-542d-48e2-8bc0-b63f98c56a4a
da_vector([1.,1.,1.,1.,1.,1.])

# ╔═╡ 2a08d783-6670-498a-b79b-d0e380c2c4e4


# ╔═╡ ec52b24b-53ea-4aba-841d-e73bed07de3d


# ╔═╡ 347a6102-01db-4779-8a36-63083147934f


# ╔═╡ 6f17b393-4645-462a-a8da-3dfc13a810d1


# ╔═╡ 8adb4e33-7470-4996-a1b2-9be09f45a06a


# ╔═╡ 2e9add50-a93b-4186-95c1-95c871f1751d


# ╔═╡ 4b422cf8-bedc-441b-862d-e7d89319bdb7


# ╔═╡ 8b8cf6b2-70ef-4c0c-b729-551637947f0c
md"The below block of code uses the vectorized derivatives above to output functions which take higher time derivatives of $\vec{\theta}$ and $I_{ij}$ as a function of $\theta$. "

# ╔═╡ 9878e729-0b7c-4684-a0fd-e0cf73465199
#=all necessary higher derivatives to calculate up to dddI=#
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

# ╔═╡ e434ffd6-ce0c-4765-8aef-564180e51309


# ╔═╡ a602a9fd-f5af-44f9-bef5-60974e100a9a



# ╔═╡ 0e80bd6c-80af-45b3-8637-22fa8bfe8a5c
#=bundling derivatives together which will prove convenient later=#
begin
    first_derivatives=[da_vector, dM_WD_vector, dM_NS_vector, dR_WD_vector, dJ_vector, dθ_vector]
    second_derivatives=[dda_vector, ddM_WD_vector, ddM_NS_vector, ddR_WD_vector, ddJ_vector, ddθ_vector]
    third_derivatives=[ddda_vector, dddM_WD_vector, dddM_NS_vector, dddR_WD_vector, dddJ_vector, dddθ_vector]
end;

# ╔═╡ a22c8171-4ca3-47e3-99c7-5e6e822983d8


# ╔═╡ 7f728824-950f-4631-bb6b-2b4654f6dc7f
#=Constructs a symmetric 2x2 tensor from a plus and times polarization=#
function tensor_plus_times(plus::Float64, times::Float64)::Matrix{Float64}
    return([plus times; times -plus        
    ])
end

# ╔═╡ 2aaf2282-5dd7-4f93-9e32-b7a45b5dc3ff


# ╔═╡ 50b7dfa1-0261-4de7-89a6-4bf1d6390310


# ╔═╡ 5098b622-aacf-4df1-8341-637e8eac60e6
md"The below block of code is central to the interface between the mass transfer and radiation reaction effects. It also sits in the innermost loop of our integrator. This is a good place to peer review."

# ╔═╡ 6208496a-07a7-4f27-9994-98a3347f018b
#=Central Engine of the integrator=#
begin
	"""Calculates the second derivative of the quadrupole moment `t_ddI`  due to kinematics of keplerian motion and mass transfer `dtheta`. Also returns higher derivatives of theta useful for calculating the third derivative of the quadrupole moment."""
	function ddI_helper(theta::AbstractVector{T}) where {T}
	    #calulates a list of values which will be helpful in the next two functions
	    dtheta=zeros(6)
	    #calculating first time derivatives due to MT only
	    for i in 1:6
	        dtheta[i]=first_derivatives[i](theta)
	    end

		theta_dtheta=vcat(theta,dtheta)
	    ddtheta=zeros(6)
	    #calculating second time derivatives due to MT only
	    for i in 1:6
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
	    for i in 1:6
	        dddtheta[i]=third_derivatives[i](theta_ddtheta)
	    end

		
	    
	    theta_dddtheta=vcat(theta_ddtheta,dddtheta)
	
	    t_dddI=tensor_plus_times(dddI_xx_vector(theta_dddtheta),
	        dddI_xy_vector(theta_dddtheta)
	    )
	    return(dtheta,t_ddI,t_dddI)
	end
end;

# ╔═╡ 8482ce78-48b8-4514-ba28-5d2ecf4f0e4e


# ╔═╡ 850fa701-4427-4f4d-bd8d-369135841d07


# ╔═╡ f9aaaa54-1e7a-46e6-954e-a6e20f637931


# ╔═╡ 2d9276d6-cf9c-439e-bd66-9220fea56dde


# ╔═╡ 2c6982e2-ee7d-4d2b-8f59-0285d85f1063


# ╔═╡ c9a96a65-a4f8-4440-840e-5467793dad63


# ╔═╡ 6a5b3c1f-f8b6-4578-915f-e35be805e608


# ╔═╡ 987850ef-c086-4b16-8401-d11c753b800d


# ╔═╡ bc193a3e-052e-4860-854d-398fa08deeec


# ╔═╡ 120ff400-43a9-4be6-b1a8-2017ae0b111a
""" An "in place" implementation of the full time derivative of theta including effects from both mass transfer and radiation reaction.""" 
function acceleration!(du::Vector{Float64}, theta::Vector{Float64},p,t)
    dtheta,t_ddI,t_dddI=ddI_dddI_from_theta(theta)

    #forces due to RR
    t_dJ_rr=dJ_rr(t_ddI,t_dddI)
    t_da_rr=2.0 * theta[1] * t_dJ_rr / theta[5]

    #adding RR forces to MT forces
    dtheta[1]+=t_da_rr
    dtheta[5]+=t_da_rr

    #mutate du
    du .= dtheta
	#println(du)
    #return(dtheta)
end;

# ╔═╡ f6882fb5-f97e-4898-855e-23b8ccc0f205
csv

# ╔═╡ 2f7ca955-4d59-483e-a86a-8363fb5795d6


# ╔═╡ f4650ecc-65c1-43d3-8d6a-676d3677df38


# ╔═╡ 947f2460-c3e0-444e-8a8c-4a3bdea80793


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
	function initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2.)::Vector{Float64}
	    #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
	    #M_NS=1.35 M_\odot is typical NS mass
	    @assert M_WD<1.4 #Chandrasekar mass, beyond which we cannot calculate R_WD(see initial function)
	    a=RL_contact(M_WD,M_WD,mu_e=mu_e)
	    J=circular_J(a,M_WD,M_WD)
	    R_WD=R0_WD(M_WD)
	    return([a, M_WD, M_NS, R_WD, J, 0.])
	end
end

# ╔═╡ cbf1359d-0d47-4f3b-a712-7ed8308b19ff
initial_circular_orbit(1.)

# ╔═╡ d22de5fc-4b5e-4a60-9e27-9879ba7f1f1f


# ╔═╡ e33d1f4e-7e39-456b-88cb-e70f3bb4dec9


# ╔═╡ 84227259-c32e-431d-8a9c-9bc1e0af252b


# ╔═╡ 1b671c3f-3fad-40aa-bf2e-e1768b6838ca


# ╔═╡ 18878cef-4449-4024-b3f3-acc108c11714


# ╔═╡ c599465d-fbef-4db4-a0ca-07edc6833525


# ╔═╡ 8032f016-4aca-4955-ab1f-2bb6c84c5713


# ╔═╡ 0b1e0b8e-b863-431f-b857-40ea9dbe37fd


# ╔═╡ 7e3f17c1-7b35-4a84-b8a5-36817549bda4


# ╔═╡ a7a7b50c-d555-4a1a-9197-963fac13aca5
md"Below is the heart of the integration process"

# ╔═╡ a8f954bf-b8d6-45f1-b836-d9a4c977f7db

begin
	"""Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken, the updated theta parameters, and the instantaneous second derivative of the quadupole moment at the end of the timestep."""
	function one_step(theta::Vector{Float64}; dt_period::Float64=1e-2, dt::Float64=-1.)::Tuple{Float64, Vector{Float64}, Array{Float64,2}}
		@assert length(theta)==6
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

	"""Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. """
	function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
	        orbit_limit::Float64=1e5, step_limit::Int=10^7)::AbstractMatrix{<:Real}
	    #returns time series second derivatives of the quadrupole moment
	    t=0.
	    theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
	    #ts=zeros(step_limit)
	    #ddIs=zeros(step_limit,2,2)
		ts_plus_times=zeros(step_limit,3)
		#print(ts_plus_times)
	    i=0
	    
	    while theta[6]<2*π*orbit_limit && i<step_limit && theta[2]>0
	        #integrate until time or step limit are elapsed or WD completely dissipates
			#println(theta)
			#println(isa(theta,Vector{Float64}))
	        dt,t_theta,t_ddI=one_step(theta; dt_period=1e-2)
	        i+=1
	        t+=dt
			#println(t," ", dt, " ", t_theta[6])
	        

			t_ddI=ddI_from_theta(theta)
			#print(t_ddI)

			#print([t, t_ddI[1,1], t_ddI[1,2]])

			#print(ts_plus_times[i])
			ts_plus_times[i,:].=[t, t_ddI[1,1], t_ddI[1,2] ]
	
	        theta=t_theta
	        #ddIs[i, :, :].=ddI_from_theta(theta)

	    end
		
	    return(ts_plus_times)
	end
end;

# ╔═╡ 84a15e51-0078-40b8-ad18-693d438eb70f


# ╔═╡ 7d40fd08-e841-401e-908f-c765d3a01f78


# ╔═╡ 56fd8555-eddb-4ad1-8ea7-273a12c49867


# ╔═╡ 3b8d5316-6342-40af-b636-d0a837376686


# ╔═╡ 2413dd9c-bd64-4923-a5c0-8471267401b8


# ╔═╡ 3aafb99a-6a64-437c-b6d1-3c7e2d63746b


# ╔═╡ c3affd3d-7d98-445d-9f8f-ae33a264dda3


# ╔═╡ e748b1fb-36a5-452c-bfdf-cc20c328575d


# ╔═╡ 63a48513-978d-47ce-afca-d58b4b9c3b2b


# ╔═╡ 07e981af-ab4d-40ff-a47e-3616e24c4596


# ╔═╡ 126d40bf-0cc2-4c7a-b07a-46bad76b011c
begin
	test_zeros=zeros(10,3)
	println(test_zeros)
	println(test_zeros[1,:])
	test_zeros[1,:].=[1,2,3]
end	

# ╔═╡ 5cea1cab-aab4-438f-a285-f4468081cb25


# ╔═╡ d7c80b56-73e2-4741-9381-23b123b7bcf0


# ╔═╡ 29342c7d-7fb5-487c-8dfc-e843d07a78e3


# ╔═╡ f58e96c8-1b26-410f-aa8c-09c26dbe2517


# ╔═╡ e6f10232-59fb-4cb4-89d2-2134b0e96a60


# ╔═╡ 1ccba8c4-964b-40c5-81e1-273d8b2afe70


# ╔═╡ dfeb7a2f-cb87-43a7-9d9f-b068d1ae509f


# ╔═╡ 094ae836-f7fd-4423-901e-89479ee28760


# ╔═╡ 99d6e037-b31a-4d47-a861-d6626128845d


# ╔═╡ 411fe157-b5f6-4961-b9b2-aacaa773253f


# ╔═╡ ea2cc713-c5e8-4d3a-bb8d-1bb868fb67ad


# ╔═╡ f00501ef-b2f0-4ca6-8225-faa31292f12f


# ╔═╡ 858944c4-387d-4880-b22a-c24561f1305c


# ╔═╡ d0bbabe3-b785-4559-b2a3-65a181de6448


# ╔═╡ 907e110a-4c66-471b-b3dc-f38f2b9d3e63


# ╔═╡ 680b68ef-660e-4e11-8129-6a1d16f24476


# ╔═╡ 7d3e0971-ebf6-4897-bbc0-c3ce5befd187


# ╔═╡ 8a28524d-7e70-40ec-aafb-53ce211fb052


# ╔═╡ 91fbe086-cf7f-433a-8df6-afbd8eebf6cc
initial_circular_orbit(1.0)

# ╔═╡ c407c3e2-ef68-4a4c-8232-029b56e44aff
period(2.35,14088.7)

# ╔═╡ bb9acbb3-9397-4d27-8bbc-21faaf415ea3


# ╔═╡ a49ee332-de38-4fd3-8d9b-04a8134f5a13
begin
	one_step([14088.687124539354, 1.0, 1.35, 5338.492628530513, 4.323851039786607e7, 0.0]; dt_period=1e-2, dt=-1.)
	println("AAAAHHHH")
end;

# ╔═╡ 828ad2f8-7d1d-45b8-9263-67ccb673cee0


# ╔═╡ 24a7185e-5dcd-47d2-b299-85bb5b210bd4


# ╔═╡ 40c0228b-c30b-42a4-b3c6-3b10d5258da8


# ╔═╡ 33ebb6b8-ba01-4941-bc51-45dec91e8bf2


# ╔═╡ 0971bd47-9d39-4029-9d39-68be71109c66


# ╔═╡ 6fdeb509-7f2a-4294-a567-05e09113b5e2


# ╔═╡ 48a16950-79e2-4bcf-9d33-83e85689e8ee


# ╔═╡ 9864a599-fb19-41a9-af26-5a3a0050975d


# ╔═╡ 5f820ec7-139e-4cb9-ab76-3be6fe0c438f


# ╔═╡ 172842d6-e181-4f86-8a2d-c871af8582fb


# ╔═╡ 1edf5071-1d21-4101-983e-3effbf694292


# ╔═╡ f9338a40-2677-4137-a2cf-da327dd157db


# ╔═╡ 1e16edf8-c038-461e-bc29-f8ed81c7caca


# ╔═╡ 91e4fc89-1284-4e31-8c97-93121782ea03
md"Below is an integration check--success depends on each function from victor.jl and tristan.jl being compatible with each other. We are currently failing this test because we are not successfully integrating victor.jl and tristan.jl (see the final line in the first block of code in this notebook."

# ╔═╡ 8f914e37-11f2-4620-b469-abc005e98312
begin
	
	test_integrate = integrate(1.0,step_limit=10)
	

	@test isa(test_integrate, AbstractMatrix{Float64})
end;

# ╔═╡ 8ad974b5-ca88-41e5-813e-ebfe86d50197


# ╔═╡ 86caa09c-5683-407b-83e8-786ba4e50116


# ╔═╡ a0869e2e-22b4-422c-904d-a57277467925
test_integrate

# ╔═╡ aff9a7ca-d8cd-4866-86f6-49d0397b0fff
begin
	test_df=DataFrame(test_integrate, :auto)
	CSV.write("Outputs/test_csv.csv", test_df)
end

# ╔═╡ adf1800e-3fea-4f8a-8255-44248a041ef4
test_read=CSV.read("Outputs/test_csv.csv", DataFrame)

# ╔═╡ f8a70c92-3d1f-48a0-a146-7f0b02e76cff
test_read

# ╔═╡ 67ee3227-0143-4575-bad0-0dbd73b812b7


# ╔═╡ b7ba9bd8-d96a-4f87-8b8b-0ec46c60fa10


# ╔═╡ c0c71b33-06aa-41ba-8a5f-ad9026fa7daa


# ╔═╡ 5fd7a5da-3505-4057-8947-7e2a061d147c
test_integrate[1]

# ╔═╡ 9526b140-ea6b-4feb-967d-16cc44d0a3af


# ╔═╡ 77d5f408-ca5d-46af-a092-7834e49a3a8c


# ╔═╡ 2ee86e92-bfd3-4e8a-b1c6-e5fb563dc32b
test_integrate[2]

# ╔═╡ ec8c7cd5-90c7-4465-ba27-f2eaedcdfc85


# ╔═╡ 3f208986-9ef2-4a0e-a89c-ce13c46a7183


# ╔═╡ 18d49c62-f3d1-4d0c-b095-fdef4a319d2d


# ╔═╡ 27eabe5a-bc40-4314-a7bd-e6993657ed1b


# ╔═╡ cd2657d4-adc8-4732-8abb-c6148f8211c6


# ╔═╡ 3ad72f7d-4c5e-416f-bd0f-bf1b04dc05b8


# ╔═╡ a1faafde-cafa-4392-bb61-e1a8af311b28


# ╔═╡ 8fbfc816-261c-4585-9f60-3740a9c68149


# ╔═╡ af8baf06-c7ce-4a17-8443-f7d988fb6257


# ╔═╡ fd2f6fb9-031c-44a4-b230-e36088190b17


# ╔═╡ 9fff71f3-9c7d-4060-a50a-cd3ec475f27c


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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
PlutoLinks = "0ff47ea0-7a50-410d-8455-4348d5de0420"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
ProfileCanvas = "efd6af41-a80b-495e-886c-e51b0c7d77a3"
SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
BenchmarkTools = "~1.6.0"
CSV = "~0.10.15"
DataFrames = "~1.8.0"
ForwardDiff = "~1.2.1"
OrdinaryDiffEq = "~6.102.1"
PlutoLinks = "~0.1.6"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.71"
ProfileCanvas = "~0.1.6"
SciMLBase = "~2.120.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "628a2aae8bddd1947218aff5ef12631559b21207"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "dbd8c3bbbdbb5c2778f85f4422c39960eac65a42"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.20.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "355ab2d61069927d4247cd69ad0e1f140b31e30d"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.12.0"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "e38fbc49a620f5d0b660d7f543db1009fe0f8336"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.6.0"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "57f9f59bec88ce80cd3e46efc39f126395789bb0"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "980f01d6d3283b3dbdfd7ed89405f96b7256ad57"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "c967271c27a95160e30432e011b58f42cd7501b5"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "087632db966c90079a5534e4147afea9136ca39a"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.190.2"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseMooncakeExt = "Mooncake"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "cee1700673af54db57bd1c7fb834ad4ff31309a0"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.8"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EnzymeCore]]
git-tree-sha1 = "e059db5d02720ae826445f5ce2fdfb3d53236b87"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.14"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cae251c76f353e32d32d76fae2fea655eab652af"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.27.0"
weakdeps = ["StaticArrays"]

    [deps.ExponentialUtilities.extensions]
    ExponentialUtilitiesStaticArraysExt = "StaticArrays"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "fd923962364b645f3719855c88f7074413a6ad92"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.0.2"

[[deps.FastPower]]
git-tree-sha1 = "5f7afd4b1a3969dc34d692da2ed856047325b06e"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.3"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerMooncakeExt = "Mooncake"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "31fd32af86234b6b71add76229d53129aa1b87a9"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.28.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "dc41303865a16274ecb8450c220021ce1e0cf05f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.2.1"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "a694e2a57394e409f7a11ee0977362a9fafcb8c7"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.6"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "d8337622fe53c05d16f031df24daf0270e53bc64"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.5"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "d1fc961038207e43982851e57ee257adc37be5e8"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.2"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "79ee64f6ba0a5a49930f51c86f60d7526b5e12c8"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.8.0"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "c4167314fda6da28522599b904bcd3f44709161f"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.42.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "e24491cb83551e44a69b9106c50666dea9d953ab"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.4"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "1d091cfece012662b06d25c792b3a43a0804c47b"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.12.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI", "SparseArrays"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "9dba8e7ccfaf4c10b3a3b4cc52abf639f8558efd"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "2.0.0"

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseChainRulesCoreExt = "ChainRulesCore"
    NonlinearSolveBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseMooncakeExt = "Mooncake"
    NonlinearSolveBaseReverseDiffExt = "ReverseDiff"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"
    NonlinearSolveBaseTrackerExt = "Tracker"

    [deps.NonlinearSolveBase.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LineSearch = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
    LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "01c48c37ba47721ec6489a1668ab3bb1f5b603c0"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.9.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "a233b7bbb32426170b238e2e2b82b0f9f1c5caba"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.10.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "139bf9211930a829703481b3ffd86b1ab309cd07"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.5.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "89cd4e81d7a668f8858fba6779212f41a0360260"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.102.1"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "09aae1486c767caa6bce9de892455cbdf5a6fbc8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.5.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "ce8db53fd1e4e41c020fd53961e7314f75e4c21c"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.10.1"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "688c717e8eee84dcfef02ddf71c8028215b329ca"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.34.0"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport", "SciMLBase"]
git-tree-sha1 = "7d5ddeee97e1bdcc848f1397cbc3d03bd57f33e7"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.8.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "320b5f3e4e61ca0ad863c63c803f69973ba6efce"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.16.1"
weakdeps = ["SparseArrays"]

    [deps.OrdinaryDiffEqDifferentiation.extensions]
    OrdinaryDiffEqDifferentiationSparseArraysExt = "SparseArrays"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "4c0633f587395d7aaec0679dc649eb03fcc74e73"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.4.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "3b81416ff11e55ea0ae7b449efc818256d9d450b"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.8.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "9e1b11cf448a2c1bca640103c1c848a20aa2f967"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.9.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "b968d66de3de5ffcf18544bc202ca792bad20710"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.16.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "815b54211201ec42b8829e0275ab3c9632d16cbe"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.4.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "fe750e4b8c1b1b9e1c1319ff2e052e83ad57b3ac"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.5.0"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "42096f72136078fa02804515f1748ddeb1f0d47d"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.5.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a5dcd75959dada0005b1707a5ca9359faa1734ba"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.7.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "925fc0136e8128fd19abf126e9358ec1f997390f"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "3cc4987c8e4725276b55a52e08b56ded4862917e"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "e6bd0a7fb6643a57b06a90415608a81aaf7bd772"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.7.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "b05226afc8fa6b8fc6f2258a89987b4f5bd0db4e"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.14.1"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "c90aa7fa0d725472c4098096adf6a08266c2f682"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.4.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "9d599d2eafdf74ab26ea6bf3feb28183a2ade143"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.6.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport", "SciMLBase"]
git-tree-sha1 = "8e35132689133255be6d63df4190b5fc97b6cf2b"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.4.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "63fb643a956b27cd0e33a3c6d910c3c118082e0f"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.4.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a31c41f9dbea7c7179c6e544c25c7e144d63868c"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.5.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "f34bc2f58656843596d09a4c4de8c20724ebc2f1"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.18.1"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "20caa72c004414435fb5769fadb711e96ed5bcd4"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.7.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "3bce87977264916bd92455754ab336faec68bf8a"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.7.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqStabilizedRK", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "75abe7462f4b0b2a2463bb512c8a5458bbd39185"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.6.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "7e94d3d1b3528b4bcf9e0248198ee0a2fd65a697"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.4.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "e8dd5ab225287947016dc144a5ded1fb83885638"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.7.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "778c7d379265f17f40dbe9aaa6f6a2a08bc7fa3e"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.5.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "185578fa7c38119d4318326f9375f1cba0f0ce53"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.6.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "c05b4c6325262152483a1ecb6c69846d2e01727b"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.34"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "60f23271568238b236899d00fc290e43d1f2658f"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.0.10"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.ProfileCanvas]]
deps = ["Base64", "JSON", "Pkg", "Profile", "REPL"]
git-tree-sha1 = "e42571ce9a614c2fbebcaa8aab23bbf8865c624e"
uuid = "efd6af41-a80b-495e-886c-e51b0c7d77a3"
version = "0.1.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "96bef5b9ac123fff1b379acf0303cf914aaabdfd"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.37.1"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "d852eba0cc08181083a58d5eb9dccaec3129cb03"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.9.0"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "16fa030fb4bd4df373a677eca0460c3eee791ab2"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.120.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "a273b291c90909ba6fe08402dd68e09aae423008"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.11"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "024d829102878141aaee5cf8f8288bcabd2f57a0"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.7.2"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "8825064775bf4ae0f22d04ea63979d8c868fd510"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.9.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "9de43e0b9b976f1019bf7a879a686c4514520078"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.21"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "1e44e7b1dbb5249876d84c32466f8988a6b41bbb"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

    [deps.StaticArrayInterface.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "83151ba8065a73f53ca2ae98bc7274d817aa30f2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.8"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "617400a198bd433f921ca2a4e89999f835dd3fde"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.45"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═6c59b881-5668-4894-950b-2c7bee5367ef
# ╠═511606a0-8594-45a3-af99-29a3ffbebe68
# ╠═8edaf48b-d27d-43c9-8a0c-e149654b6cd0
# ╠═0cac30f1-9101-4ba3-accc-52621bc1d16f
# ╠═c79771dc-4bf2-4c03-9f84-dd6532c98bfd
# ╠═33f9cc9a-a207-4a09-bfa3-0f54546c05c3
# ╠═e2833dab-01a3-4d20-a42d-7d182f9e610d
# ╠═1db7cf48-6630-4701-b163-c2fcc632a050
# ╠═5519a048-4a4a-4091-af9e-8b2dc444775b
# ╠═c13a77be-813e-413f-b0b2-23cf7b5f309a
# ╠═14fedf73-b9ca-4a37-a86c-e6f4c0a59ece
# ╠═abca47ab-091a-4677-8d5c-b13572bb7baa
# ╠═b97703e5-e4f3-46a9-8819-81568e8582d6
# ╠═4bef6030-10c5-4497-9ec2-edabee8432d1
# ╠═17be5889-9440-4adb-965f-17de4875a922
# ╠═ce80bb70-d60e-4bdf-8f7c-16421fd90138
# ╠═fa5a36f9-b1e2-4461-bb18-a3172f4b7af7
# ╠═421df964-fd81-471d-bff5-570f3f409f53
# ╠═2ea8a9da-f653-41f8-97df-f2909f7ce1a2
# ╠═25d4f52f-7f37-41c9-9230-37dfb0e71f22
# ╠═104e6849-3bae-467a-8630-57297b116314
# ╠═3ef782fd-4d6a-4c66-a976-5e995e2d703d
# ╠═13bef496-5bb6-4dbe-a6e4-d18e664dfc55
# ╠═9a2f1c1e-46ca-41b8-9d11-d165bb86a918
# ╠═3e28aef1-3662-47b6-8d59-9a1c97313d5b
# ╠═190a4e03-541a-4857-879d-7093afdfddc2
# ╠═d3985c1e-dde5-4700-b948-b6fd04c80819
# ╠═40d859c3-a08e-49a2-8b8b-34ad7a0a6fa3
# ╠═35d54cfa-7d0a-4a48-94bc-a45327b0abcc
# ╠═1620f45c-a715-432c-8a2d-f673bb022e2c
# ╠═96c169cc-f4d7-4885-ab1f-a692f5f39f91
# ╠═d0cc1fbf-6828-4bf5-8bb6-1f9d97161a24
# ╠═1ffcf92d-958a-44a2-88cc-1c7ab262925b
# ╠═8aca584e-0052-4e4d-9e30-823204f9595e
# ╠═fc8da2c6-e7b8-47cd-af05-b13e74909dda
# ╠═e4452b17-0722-4ea0-804b-7f2797e17a82
# ╠═ac043976-9ee2-4633-915b-1043a4771164
# ╠═647e477c-681f-43e7-9794-38d19bca6098
# ╠═d06e9c81-4436-4054-a7cc-c85b7ee13df9
# ╠═07ca0c95-f9b3-42ad-8a49-3d305002b20f
# ╠═4f1118f3-e8f8-479f-819a-db201f36e447
# ╠═86420e50-5f10-43aa-ac65-a8708be68dc6
# ╠═6d4270fb-21e4-4e8c-b7e4-0e7bc83870c6
# ╠═94098f26-fc8d-4df5-84fe-9ebb7b8c994a
# ╠═256492bf-b826-49aa-b051-55991228703c
# ╠═a3e1afb5-0a7b-41aa-b50f-dc959a36a1e6
# ╠═69081296-10a7-4628-a62d-b116548bdbae
# ╠═bfbaf692-db1f-41f9-a147-95ccb8de9e1e
# ╠═1de07299-8623-49fe-bba3-136dc45e7c82
# ╠═ea61a15f-1982-4469-8e0b-22f309d1bef5
# ╠═89014777-e6e0-406a-863c-72a310e37f7e
# ╠═45d20cd8-a5e6-4709-b309-bfb3da22ed67
# ╠═b2dcf719-078c-4ec0-bde6-243972cfc957
# ╠═4135feba-248d-4880-9e05-e3e85dae9963
# ╠═2fb18e7f-30eb-4dff-b30f-f5086ec3e02d
# ╠═f784cf7a-6cf3-48b0-b705-acf708b75db8
# ╠═d43a93a0-3f71-4d7f-b548-6057a1513f87
# ╠═d45576f6-f2ad-47e5-a4a7-6bf6f95d42c9
# ╠═557daf12-69f8-4f4f-894f-6bc435fa91ab
# ╠═df23f515-ffc2-4965-b695-c3875c5bf29d
# ╠═760d0fda-81d0-4d36-9e3d-46912d0ddfa4
# ╠═6ac10460-5123-487d-a55d-6183de743a5e
# ╠═d93a9205-0f8c-4b47-af4e-5ad0b6437353
# ╠═344fdbd0-96e3-4a34-8d9f-3b6c34cee2b3
# ╠═cbf1359d-0d47-4f3b-a712-7ed8308b19ff
# ╠═5512bb36-6de6-41e6-93a4-d736f8b69b3d
# ╠═a1e4b8ec-de2e-42ab-8d52-2f62cf5d92dd
# ╠═c16514a9-c297-451b-9e6d-1d21d15d0f28
# ╠═d9e0ca74-8f4b-4d10-8255-7a1dae0d49d7
# ╠═7fc2e8cf-556f-421f-a993-34e845c4ed3f
# ╠═0757c720-c651-44bc-aca6-20798890e01b
# ╠═5a5761f7-eb8e-419f-8bb6-15780de928d2
# ╠═8f6900c9-352f-4deb-82e2-1bd99faffda4
# ╠═bb73d75d-eeeb-4b38-8db2-5f9cd4a440c3
# ╠═15095435-54c0-4d8d-91c0-c278edd57700
# ╠═870a48df-44ba-4585-b422-c92a84737a12
# ╠═2e25bf2a-d4c0-41c6-b4f3-1fff15bd6032
# ╠═7ffbc67d-2ad6-4ca2-8e52-7bcd81089da1
# ╠═a18ec6f3-6036-4e0e-b8f2-6156e0497e8c
# ╠═69640bb5-07a4-4257-9eea-6f16dfcce74c
# ╠═7faf665d-0bdc-499d-8cc6-72ebd7da30b7
# ╠═c42a4b21-75f4-4246-a34c-c51fdd752d0a
# ╠═94e9ebac-2ffa-4e76-be6c-273f402064ae
# ╠═57900831-e4f4-4557-bacc-422af3cc2d2b
# ╠═02584e20-271e-4263-a0b9-33889bc18ae6
# ╠═dd2c1373-c195-49b1-8f00-18894900bd85
# ╠═c69bd40e-49a0-4f5c-b516-6d5f774dd75d
# ╠═2efd9365-9bcf-41ea-9e77-9346bf498bde
# ╠═6b59b741-3cfa-4db2-afce-e6f23ab1a06a
# ╠═4bbe25ad-bce2-4305-bd83-6dba2dc1f28c
# ╠═e90a2a13-e9af-4160-a307-cc7df36e280e
# ╠═aad547b0-ff39-4c6c-bd40-d47d132e3f4d
# ╠═d877e869-f0ec-4b02-9bdc-8f8d6d4e1c6c
# ╠═deb8fbb9-e4f9-44e3-87cf-e226319b8a2e
# ╠═3fb94d9d-948f-4620-8621-e0f3d11e99cf
# ╠═84c2d206-599a-42e1-ba61-6bc06cb62f65
# ╠═7f29652e-29d7-43dc-8ad3-a1dae2f39991
# ╠═69e12b46-10d9-497b-9d2c-7b54e1658d7f
# ╠═61544091-4d40-434b-bf52-de66d8c35c7a
# ╠═25998cc0-96e9-476c-8853-d102a58c0334
# ╠═0a57469d-418f-47b0-bf56-828d4bc529e1
# ╠═694fe531-e34e-4c27-94cb-ae3f1e14aa90
# ╠═80ee41e0-5978-47fe-8554-16477fb317d6
# ╠═ba3f4fbe-1eb3-4040-8616-28f926bb6243
# ╠═fdd55e22-89d2-465a-80b8-5c59f5f7129b
# ╠═92a94cfd-a3ff-4d62-af3b-dc8b0b57b0ba
# ╠═28bdbf35-455c-4e01-a1bb-05c862abe361
# ╠═72720c31-b465-4755-b070-9f34311a0d97
# ╠═e1591006-147d-484b-8559-763e239af504
# ╠═57002625-8ba1-43a1-b2e3-0094fa53c841
# ╠═5e8ec2f8-d73a-4a81-86dc-bb9c3b538caa
# ╠═c1636ad3-16b4-43d1-b2bb-993885d11818
# ╠═44ba797c-48c8-45e2-87af-a6bf102066a3
# ╠═951f256e-c418-4f39-b4aa-0f4080caae8b
# ╠═a238a819-a6ef-479c-912a-2ddc5812ce66
# ╠═72dd65d6-0b35-4d2e-a359-3abe9aa39c5f
# ╠═4606dfbb-e601-4da5-8ed0-e7cf1e5be3db
# ╠═bd182ff6-48cc-445e-b067-e387b88ff08b
# ╠═4488d59c-5ef0-4d3e-8695-161bfefaf51c
# ╠═6951b817-887e-46b7-9e5c-e4840ffd6cf4
# ╠═d71448b4-705e-468a-91b5-d3c87b599316
# ╠═a3315e6b-846c-4f8b-96bf-aa2519d44ae7
# ╠═4a299a1e-afeb-4c7c-82a9-a3140a093719
# ╠═1efcd158-2476-43ca-9474-66ec7e391aa0
# ╠═f52559d9-35d7-4464-bd88-26419fb5354f
# ╠═d95b1547-717e-4ba4-8a1c-edea66926550
# ╠═11d2af20-64fd-41de-b33c-76122c047bcb
# ╠═621d6981-fed6-4665-85a2-f5df3a06519e
# ╠═05576550-3ffb-4d2e-b89e-e4d00ed54d71
# ╠═509ec580-b9b1-42dd-ab81-b72d1e1f7293
# ╠═19be5fd8-ea30-4cb7-93a4-cc2171651194
# ╠═5bdbc460-4df3-4d64-87d2-7820997d2317
# ╠═b889ebe5-e259-4208-8685-cc32df37f58c
# ╠═95d0ce5f-fe9a-4cd2-bb0f-90eb4a5e1f39
# ╠═d41d2eb7-701e-4893-8116-7d477effeef4
# ╠═de0d4bc7-e2b8-4812-8b28-43bda221b17f
# ╠═31736097-efce-4c31-85ab-d657cb97822a
# ╠═19a2498d-e8c1-4394-bc7c-9dd2e6033b16
# ╠═96a595df-5e71-4d63-800f-f600ef4b591a
# ╠═e681a7b9-9996-4105-ac4c-f1bd880f8235
# ╠═7d4b6bb6-2fcc-4669-9c4d-369b10b2644b
# ╠═7fcacf8f-b82f-48f0-87ab-9f11b7ac936f
# ╠═f9e6f40b-ebb8-4711-99b2-6852b6a9c223
# ╠═ec453c88-b887-472e-a98c-f96c3b896f08
# ╠═9d610030-0721-440d-a7dd-ea6007249228
# ╠═24cb207f-cd3a-4e50-ad02-6318c607738c
# ╠═731205f1-64f9-4dd5-809f-a1cd20eabd47
# ╠═528d99bb-0800-41d7-812b-8685a18c2ec0
# ╠═3b211714-36b0-4b0d-9b62-6ae2939f66d9
# ╠═bf13e120-01a9-4026-90ee-29576bfa784c
# ╠═680002e1-a8e8-4fec-a37c-e6c92d8a6d2e
# ╠═59b217a8-0b07-442e-bcd4-7fda42d5d03a
# ╠═d0db0f15-0632-4f3b-bdc0-af1342ae03a7
# ╠═2f9e1da8-3522-4aad-af07-00d13602b610
# ╠═fd1ffa20-5780-4ee5-8f6d-248c02a1656e
# ╠═736914e5-11e5-40bb-8792-236e132699d5
# ╠═47c8e194-843f-4648-95c7-ad08a96cba19
# ╠═42d106c9-1b6f-42ed-966b-4c5fb40b7f5f
# ╠═2f6194b0-615f-4342-a49a-7364756b4d01
# ╠═cf21fed8-294e-40ee-b29a-b1e0b673b834
# ╠═a00e5b5a-488e-413e-aea1-f404dec591b5
# ╠═4ba2c0d2-c006-4a2d-ae05-48fc40948cf6
# ╠═e2064499-24de-4c0a-a975-289dcf61b361
# ╠═ecce3a10-b519-49d9-ab59-bd33926b5e3d
# ╠═ad39dbb0-a043-4ed9-ada1-5a59612f57ec
# ╠═e70e5b74-542d-48e2-8bc0-b63f98c56a4a
# ╠═2a08d783-6670-498a-b79b-d0e380c2c4e4
# ╠═ec52b24b-53ea-4aba-841d-e73bed07de3d
# ╠═347a6102-01db-4779-8a36-63083147934f
# ╠═6f17b393-4645-462a-a8da-3dfc13a810d1
# ╠═8adb4e33-7470-4996-a1b2-9be09f45a06a
# ╠═2e9add50-a93b-4186-95c1-95c871f1751d
# ╠═4b422cf8-bedc-441b-862d-e7d89319bdb7
# ╠═8b8cf6b2-70ef-4c0c-b729-551637947f0c
# ╠═9878e729-0b7c-4684-a0fd-e0cf73465199
# ╠═e434ffd6-ce0c-4765-8aef-564180e51309
# ╠═a602a9fd-f5af-44f9-bef5-60974e100a9a
# ╠═0e80bd6c-80af-45b3-8637-22fa8bfe8a5c
# ╠═a22c8171-4ca3-47e3-99c7-5e6e822983d8
# ╠═7f728824-950f-4631-bb6b-2b4654f6dc7f
# ╠═2aaf2282-5dd7-4f93-9e32-b7a45b5dc3ff
# ╠═50b7dfa1-0261-4de7-89a6-4bf1d6390310
# ╠═5098b622-aacf-4df1-8341-637e8eac60e6
# ╠═6208496a-07a7-4f27-9994-98a3347f018b
# ╠═8482ce78-48b8-4514-ba28-5d2ecf4f0e4e
# ╠═850fa701-4427-4f4d-bd8d-369135841d07
# ╠═f9aaaa54-1e7a-46e6-954e-a6e20f637931
# ╠═2d9276d6-cf9c-439e-bd66-9220fea56dde
# ╠═2c6982e2-ee7d-4d2b-8f59-0285d85f1063
# ╠═c9a96a65-a4f8-4440-840e-5467793dad63
# ╠═6a5b3c1f-f8b6-4578-915f-e35be805e608
# ╠═987850ef-c086-4b16-8401-d11c753b800d
# ╠═bc193a3e-052e-4860-854d-398fa08deeec
# ╠═120ff400-43a9-4be6-b1a8-2017ae0b111a
# ╠═f6882fb5-f97e-4898-855e-23b8ccc0f205
# ╠═2f7ca955-4d59-483e-a86a-8363fb5795d6
# ╠═f4650ecc-65c1-43d3-8d6a-676d3677df38
# ╠═947f2460-c3e0-444e-8a8c-4a3bdea80793
# ╠═52b828d3-a6b7-4aa6-abe5-5c23ef848bb2
# ╠═a63692e2-f652-4f35-b5c8-57ed131d676e
# ╠═72bfe0e8-11f2-47f6-97ba-09708af3fad3
# ╠═d22de5fc-4b5e-4a60-9e27-9879ba7f1f1f
# ╠═e33d1f4e-7e39-456b-88cb-e70f3bb4dec9
# ╠═84227259-c32e-431d-8a9c-9bc1e0af252b
# ╠═1b671c3f-3fad-40aa-bf2e-e1768b6838ca
# ╠═18878cef-4449-4024-b3f3-acc108c11714
# ╠═c599465d-fbef-4db4-a0ca-07edc6833525
# ╠═8032f016-4aca-4955-ab1f-2bb6c84c5713
# ╠═0b1e0b8e-b863-431f-b857-40ea9dbe37fd
# ╠═7e3f17c1-7b35-4a84-b8a5-36817549bda4
# ╠═a7a7b50c-d555-4a1a-9197-963fac13aca5
# ╠═a8f954bf-b8d6-45f1-b836-d9a4c977f7db
# ╠═84a15e51-0078-40b8-ad18-693d438eb70f
# ╠═7d40fd08-e841-401e-908f-c765d3a01f78
# ╠═56fd8555-eddb-4ad1-8ea7-273a12c49867
# ╠═3b8d5316-6342-40af-b636-d0a837376686
# ╠═2413dd9c-bd64-4923-a5c0-8471267401b8
# ╠═3aafb99a-6a64-437c-b6d1-3c7e2d63746b
# ╠═c3affd3d-7d98-445d-9f8f-ae33a264dda3
# ╠═e748b1fb-36a5-452c-bfdf-cc20c328575d
# ╠═63a48513-978d-47ce-afca-d58b4b9c3b2b
# ╠═07e981af-ab4d-40ff-a47e-3616e24c4596
# ╠═126d40bf-0cc2-4c7a-b07a-46bad76b011c
# ╠═5cea1cab-aab4-438f-a285-f4468081cb25
# ╠═d7c80b56-73e2-4741-9381-23b123b7bcf0
# ╠═29342c7d-7fb5-487c-8dfc-e843d07a78e3
# ╠═f58e96c8-1b26-410f-aa8c-09c26dbe2517
# ╠═e6f10232-59fb-4cb4-89d2-2134b0e96a60
# ╠═1ccba8c4-964b-40c5-81e1-273d8b2afe70
# ╠═dfeb7a2f-cb87-43a7-9d9f-b068d1ae509f
# ╠═094ae836-f7fd-4423-901e-89479ee28760
# ╠═99d6e037-b31a-4d47-a861-d6626128845d
# ╠═411fe157-b5f6-4961-b9b2-aacaa773253f
# ╠═ea2cc713-c5e8-4d3a-bb8d-1bb868fb67ad
# ╠═f00501ef-b2f0-4ca6-8225-faa31292f12f
# ╠═858944c4-387d-4880-b22a-c24561f1305c
# ╠═d0bbabe3-b785-4559-b2a3-65a181de6448
# ╠═907e110a-4c66-471b-b3dc-f38f2b9d3e63
# ╠═680b68ef-660e-4e11-8129-6a1d16f24476
# ╠═7d3e0971-ebf6-4897-bbc0-c3ce5befd187
# ╠═8a28524d-7e70-40ec-aafb-53ce211fb052
# ╠═91fbe086-cf7f-433a-8df6-afbd8eebf6cc
# ╠═c407c3e2-ef68-4a4c-8232-029b56e44aff
# ╠═bb9acbb3-9397-4d27-8bbc-21faaf415ea3
# ╠═a49ee332-de38-4fd3-8d9b-04a8134f5a13
# ╠═828ad2f8-7d1d-45b8-9263-67ccb673cee0
# ╠═24a7185e-5dcd-47d2-b299-85bb5b210bd4
# ╠═40c0228b-c30b-42a4-b3c6-3b10d5258da8
# ╠═33ebb6b8-ba01-4941-bc51-45dec91e8bf2
# ╠═0971bd47-9d39-4029-9d39-68be71109c66
# ╠═6fdeb509-7f2a-4294-a567-05e09113b5e2
# ╠═48a16950-79e2-4bcf-9d33-83e85689e8ee
# ╠═9864a599-fb19-41a9-af26-5a3a0050975d
# ╠═5f820ec7-139e-4cb9-ab76-3be6fe0c438f
# ╠═172842d6-e181-4f86-8a2d-c871af8582fb
# ╠═1edf5071-1d21-4101-983e-3effbf694292
# ╠═f9338a40-2677-4137-a2cf-da327dd157db
# ╠═1e16edf8-c038-461e-bc29-f8ed81c7caca
# ╠═91e4fc89-1284-4e31-8c97-93121782ea03
# ╠═8f914e37-11f2-4620-b469-abc005e98312
# ╠═8ad974b5-ca88-41e5-813e-ebfe86d50197
# ╠═86caa09c-5683-407b-83e8-786ba4e50116
# ╠═a0869e2e-22b4-422c-904d-a57277467925
# ╠═aff9a7ca-d8cd-4866-86f6-49d0397b0fff
# ╠═adf1800e-3fea-4f8a-8255-44248a041ef4
# ╠═f8a70c92-3d1f-48a0-a146-7f0b02e76cff
# ╠═67ee3227-0143-4575-bad0-0dbd73b812b7
# ╠═b7ba9bd8-d96a-4f87-8b8b-0ec46c60fa10
# ╠═c0c71b33-06aa-41ba-8a5f-ad9026fa7daa
# ╠═5fd7a5da-3505-4057-8947-7e2a061d147c
# ╠═9526b140-ea6b-4feb-967d-16cc44d0a3af
# ╠═77d5f408-ca5d-46af-a092-7834e49a3a8c
# ╠═2ee86e92-bfd3-4e8a-b1c6-e5fb563dc32b
# ╠═ec8c7cd5-90c7-4465-ba27-f2eaedcdfc85
# ╠═3f208986-9ef2-4a0e-a89c-ce13c46a7183
# ╠═18d49c62-f3d1-4d0c-b095-fdef4a319d2d
# ╠═27eabe5a-bc40-4314-a7bd-e6993657ed1b
# ╠═cd2657d4-adc8-4732-8abb-c6148f8211c6
# ╠═3ad72f7d-4c5e-416f-bd0f-bf1b04dc05b8
# ╠═a1faafde-cafa-4392-bb61-e1a8af311b28
# ╠═8fbfc816-261c-4585-9f60-3740a9c68149
# ╠═af8baf06-c7ce-4a17-8443-f7d988fb6257
# ╠═fd2f6fb9-031c-44a4-b230-e36088190b17
# ╠═9fff71f3-9c7d-4060-a50a-cd3ec475f27c
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
