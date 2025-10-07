using ForwardDiff
#package used for chain rules of large derivatives

#smooth approximation to the maximum function
function Boltzmann_Max(x1::Float64, x2::Float64; α::Float64=4)::Float64
    return((x1*exp(α*x1)+x2*exp(α*x2)/(exp(α*x1)+exp(α*x2))))
end
    

function dM_NS(a::Float64, R_WD::Float64, M_WD::Float64, M_NS::Float64;
        A::Float64=10, mu_e::Float64=2.0, eta_acc::Float64=0.1, α::Float64=4)::Float64
    t_dM_WD=decretion_rate(a,a,R_WD,M_WD,M_Ns,A=A)
    t_edd=eddington_rate(M_NS, mu_e=mu_e, eta_acc=eta_acc)
    return(Boltzmann_Max(-t_dM_WD,t_edd,α=α))
end

function dθ(J::Float64,a::Float64,M_WD::Float64, M_NS::Float64)::Float64
        @assert M_WD > 0 # In production-ready code, consider commenting this out
    	@assert M_NS > 0 # In production-ready code, consider commenting this out
        @assert r > 0 # In production-ready code, consider commenting this out
        return(J/(a*mu(M_WD,M_NS)))
end        
    

function I_xx(M_WD::Float64, M_NS::Float64, r::Float64, θ::Float64)::Float64
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert r > 0 # In production-ready code, consider commenting this out
    return (mu(M_WD, M_NS) * (x(r,θ))^2)
end           
            

function I_xy(M_WD::Float64, M_NS::Float64, r::Float64, θ::Float64)::Float64
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert r > 0 # In production-ready code, consider commenting this out
    return mu(M_WD, M_NS) * x(r,θ) * y(r,θ)
end                

    
#=Calculates time derivative of I_xx/I_xy from autodifferentiation=#
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


#=The following functions wrap the original derivatives in vectors to make life easier 
when taking derivatives due to mass transfer only=#    
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

function dM_R_WD_vector(theta::Vector{Float64}; mu_e=2::Float64)::Float64
    a, M_WD, M_NS, R_WD, J, θ=theta
    return( R_WD(M_WD, mu_e=mu_e))
end
    
function dJ_vector(theta::Vector{Float64}; 
               A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
               a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64)::Float64
    a, M_WD, M_NS, R_WD, J, θ=theta
    return(J_dot(J, a, a, R_WD, M_WD, M_NS, 
               A=A, mu_e=mu_e, eta_acc=eta_acc, 
               a_ring_frac=a_ring_frac, mode=mode))
end
function dθ_vector(theta::Vector{Float64})::Float64
    a, M_WD, M_NS, R_WD, J, θ=theta
    return(dθ(J,a,M_WD, M_NS))
end

function 


    

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

                    