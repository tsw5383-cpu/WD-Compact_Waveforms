using ForwardDiff
#package used for chain rules of large derivatives
using OrdinaryDiffEq
using SciMLBase
#for integrating differential equations

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


#=Same treatment, but for quadrupole moments=#
function Ixx_vector(theta::Vector{Float64})::Float64
    a, M_WD, M_NS, R_WD, J, θ=theta
    return(I_xx(M_WD, M_NS, a, θ))
end

function Ixy_vector(theta::Vector{Float64})::Float64
    a, M_WD, M_NS, R_WD, J, θ=theta
    return(I_xy(M_WD, M_NS, a, θ))
end



#=Higher derivatives=#
#=takes in values and time derivatives to calculate second derivative=#

function next_derivative(in_ftn::Function; n::Integer=1)::Function
    #returns a function which is the next higher derivative of the in_ftn
    #n is the total derivatives past the given analytic form
    function ftn(gen_theta::Vector{Float64}, gen_dtheta::Vector{Float64})::Float64
        @assert length(gen_theta)==length(gen_dtheta)==6*(1+n) 
        @assert gen_theta[7:]==gen_dtheta[:-6]
        #requires gen_theta to be [theta, dtheta, d...dtheta] for n=number of ds
        #requires gen_theta to be [dtheta, ddtheta d...ddtheta] for n+1=number of ds
        return(ftn(in_ftn, gentheta, gen_dtheta))
    end
    return ftn
end

#=all necessary higher derivatives to calculate up to dddI=#
begin
    dda=next_derivative(da_vector)
    ddM_WD=next_derivative(dM_WD_vector)
    ddM_NS=next_derivative(dM_WD_vector)
    ddR_WD=next_derivative(dR_WD_vector)
    ddJ=next_derivative(dJ_vector)
    ddθ=next_derivative(dθ_vector)

    ddda=next_derivative(dda_vector,n=2)
    dddM_WD=next_derivative(ddM_WD_vector,n=2)
    dddM_NS=next_derivative(ddM_NS_vector,n=2)
    dddR_WD=next_derivative(ddR_WD_vector,n=2)
    dddJ=next_derivative(ddJ_vector,n=2)
    dddθ=next_derivative(ddθ_vector,n=2)

    dI_xx=next_derivative(Ixx_vector)
    ddI_xx=next_derivative(Ixx_vector,n=2)
    dddI_xx=next_derivative(Ixx_vector,n=3)
    
    dI_xy=next_derivative(Ixy_vector)
    ddI_xy=next_derivative(Ixy_vector,n=2)
    dddI_xy=next_derivative(Ixy_vector,n=3)
end;    

#=bundling derivatives together which will prove convenient later=#
begin
    first_derivatives=[da_vector dM_WD_vector dM_WD_vector dR_WD_vector dJ_vector dθ_vector]
    second_derivatives=[dda_vector ddM_WD_vector ddM_WD_vector ddR_WD_vector ddJ_vector ddθ_vector]
    third_derivatives=[ddda_vector dddM_WD_vector dddM_WD_vector dddR_WD_vector dddJ_vector dddθ_vector]
end;

#=Constructs a symmetric 2x2 tensor from a plus and times polarization=#
function tensor_plus_times(plus::Float64, times::Float64)::Vector{Float64}
    return([plus times;
            times -plus        
    ])
end

#=Central Engine of the integrator=#

function ddI_dddI_from_theta(theta::Vector{Float64};ddI_only::Boolean=false)::Vector{Float64}
    #returns time derivatives of theta and second and third derivatives of 
    #the quadrupole moment due to MT effects
    
    theta_ddtheta=zeros(18)#holds theta and the first two time derivatives
    theta_ddtheta[1:6]=theta
        
    dtheta=zeros(6)
    #calculating first time derivatives due to MT only
    for i in 1:6
        dtheta[i]=first_derivatives(theta)
    end
    
    #theta_dtheta=vcat(theta,dtheta)
    theta_ddtheta[7:12]=dtheta

    ddtheta=zeros(6)
    #calculating second time derivatives due to MT only
    for i in 1:6
        ddtheta[i]=second_derivatives(theta_dddtheta[1:12])
    end

    theta_ddtheta[13:18]=ddtheta

    if not ddI_only
        dddtheta=zeros(6)
        #calculating third time derivatives due to MT only
        for i in 1:6
            dddtheta[i]=theta_derivatives(theta_dddtheta[1:18])
        end
    end


    
    #using orbital elements and first and second time derivatives 
    #to calculate second derivative of quadrupole moment
    #(for use in calculating RR acceleration)
    
    dtheta_ddtheta=vcat(dtheta,ddtheta)
    t_ddI=tensor_plus_times(ddI_xx(theta_dtheta,dtheta_ddtheta),
        ddI_xy(theta_dtheta,dtheta_ddtheta)
    )
    if ddI_only
        return(t_ddI)
    else 
        #using orbital elements and first through third time derivatives  
        #to calculate third derivative of quadrupole moment 
        #(for use in calculating RR acceleration)
        dtheta_dddtheta=vcat(dtheta_ddtheta,dddtheta)
        t_dddI=tensor_plus_times(dddI_xx(theta_ddtheta,dddtheta),
            dddI_xy(theta_ddtheta,dddtheta)
        )
        return(dtheta,t_ddI,t_dddI)
    end
end

#I am confused on how ODE solving works in Julia...
function acceleration!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Float64},t::Float64)::Vector{Float64}
    dtheta,t_ddI,t_dddI=ddI_dddI_from_theta(theta)

    #forces due to RR
    t_dJ_rr=dJ_rr(t_ddI,t_dddI)
    t_da_rr=2.0 * theta[1] * dJ_rr(I2, I3) / theta[5]

    #adding RR forces to MT forces
    dtheta[1]+=t_da_rr
    dtheta[5]+=t_da_rr
    return(dtheta)
end


#=The following functions are defined to initiate the orbit beginning at mass transfer
for initial masses, take an integration step forward calculating waveform, and integrate
until the White Dwarf dissipates completely=#

function RL_contact(M_WD::Float64,M_NS::Float64;mu_e=2::Float64)::Float64
    #returns separation at which a donor WD overflows its
    #Roche Lobe
    R1=R_WD(M_WD,mu_e=mu_e)
    RL_a1=Roche_Limit(1,M1,M2) #roche limit at separation of 1 km
    return(R1 / RL_a1)#initial separation in km
end

function circular_J(a::Float64,M_WD::Float64,M_NS::Float64)::Float64
    #returns angular momentum for a circular orbit
    return((2 * a * G * M_WD ^ 2 * M_NS ^ 2 / (M1 + M2))^(1/2))
end

function initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e=2::Float64)::Float64
    #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
    #M_NS=1.35 M_\odot is typical NS mass 
    a=initial_separation(M_WD,M_WD,frac=1e-1)
    J=circular_J(a,M_WD,M_WD)
    R1=R_WD(M_WD)
    return([a M_WD M_NS R_1 J 0])
end

#=I am confused as to whether I am treating ODEProblem correctly=#

function one_step(theta::Float64; dt_period::Float64=1e-2, dt::Float64=-1)
    if dt==-1
        dt=dt_period*period(theta[2]+theta[3], theta[1])
    end
    prob = ODEProblem(acceleration,theta,(0,dt))
    updated_theta = solve(prob, alg=KahanLi8())

    #=Calculating ddI at the end of the step=#
        
    t_ddI=ddI_dddI_from_theta(updated_theta, ddI_only=true)
    return(dt,updated_theta,t_ddI)
end

function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e=2::Float64, 
        orbit_limit::Float64=1e5, step_limit::Int=10^7)
    #returns time series second derivatives of the quadrupole moment
    t=0.
    theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
    ts=np.zeros(step_limit)
    =np.zeros(step_limit,2,2)
    i=0
    time_limit=orbit_limit*period(M_WD+M_NS,theta[1])
    
    while t<time_limit && i<step_limit && theta[2]>0
        #integrate until time or step limit are elapsed or WD completely dissipates
        dt,t_theta,t_ddI=one_step(theta)
        i+=1
        t+=dt
        ts[i]=t

        theta=t_theta
        ddIs[i, :, :].=ddI_dddI_from_theta(theta,ddI_only=true)
    end

    return(ts[1:i],ddIs[1:i, :, :])
end
    
function gravitational_waveform_ftn(i::Float64, φ::Float64, DL::Float64 )::Function
    #fixes a function to return gravitational wave given orientation and distance
    #DL is in Mpc
    @assert 0<=i<=π
    @assert 0<=φ<2π
    @assert DL<10000 #Ensure correct units for DL
    e1=[cos(i)*cos(φ)   sin(i)*sin(φ)  cos(i)  ]
    e2=[-sin(φ)  cos(φ)  0]
    es=[e1,e2]

    matrices=zeros(4,2,2)
    for i in 1:2
        for j in 1:2
            matrices[2*(i-1)+j ,:,:] = es[i] * es[j] #OUTER PRODUCT OF VECTORS
        end
    end

    factor= 2*G/c^4
    
    function ftn(ddI::Vector{Float64})::Vector{Float64}
        #given fixed outer products, calculate strain
        prod=zeros(2,2) #will end up being in transverse gauge
        for i in 1:2
            for j in 1:2
                prod[i,j] = sum(matrices[i,j]) #OUTER PRODUCT OF VECTORS    
            end
        end
        prod_TT=prod .- (prod[1,1]+prod[2,2]) .* [1 0; 0 1] #traceless
        return((factor/DL))
    
    
    


    
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

                    