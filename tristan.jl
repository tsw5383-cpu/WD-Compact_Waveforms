using ForwardDiff
#package used for chain rules of large derivatives

function dM_NS(a::Float64, R_WD::Float64, M_WD::Float64, M_NS::Float64;
        A::Float64=10, mu_e::Float64=2.0, eta_acc::Float64=0.1)
    t_dM_WD=decretion_rate(a,a,R_WD,M_WD,M_Ns,A=A)
    t_edd=eddington_rate(M_NS; mu_e=mu_e, eta_acc=eta_acc)
    return()
    

function I_xx(M_WD::Float64, M_NS::Float64, r::Float64, θ::Float64)
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert r > 0 # In production-ready code, consider commenting this out
    return mu(M_WD, M_NS) * (x(r,θ))^2

function I_xy(M_WD::Float64, M_NS::Float64, r::Float64, θ::Float64)
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert r > 0 # In production-ready code, consider commenting this out
    return mu(M_WD, M_NS) * x(r,θ) * y(r,θ)

    

function dI_xx(M_WD::Float64, M_NS::Float64, a::Float64, θ::Float64, R_WD::Float64; A::Float64=10.0)
    @assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    @assert a > 0 # In production-ready code, consider commenting this out
            
    grad=ForwardDiff.gradient(v->I_xx(v[1],v[2],v[3],v[4]), 
        [M_WD, M_NS, a, θ])#gradient of I_xx  wrt tracked parameters

    #time derivatives of  M_WD,M_NS,a,θ       
    t_dM_WD=decretion_rate(a, a, R_WD, M_WD, M_NS, A=A)
    
            
    