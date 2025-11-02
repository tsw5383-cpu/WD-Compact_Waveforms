### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

begin
	using LinearAlgebra
	using BenchmarkTools
	using Profile, ProfileCanvas 
	using PlutoUI, PlutoTeachingTools, PlutoLinks
    using Test
end;


# Constants
begin
	const G::Float64 = 1.327e11   # km^3 / M⊙ / s^2
	const c::Float64 = 2.998e5    # km/s
    const period_factor::Float64 = sqrt(4.0 * π^2 / G)
    const dJ_rr_factor::Float64 = -2.0 * G / (5.0 * c^5)
end




function R0_WD(M_WD;mu_e=2.0) 
	#@assert M_WD > 0 # In production ready code, consider commenting out assertion
	@assert mu_e > 0 # In production ready code, consider commenting out assertion
	#white dwarf radius in km as ftn of mass in M_odot and mean molecular weight/electron=2
	return(10^4 * (M_WD/0.7)^(-1/3) * (1-(M_WD/1.4)^(4/3))^(1/2) * (mu_e/2)^(-5/3))
end

function dynamical_timescale(M_WD,R_WD)
	#@assert M_WD > 0 # In production ready code, consider commenting out assertion
	#@assert R_WD > 0 # In production ready code, consider commenting out assertion
	#returns dynamical timescale for perturbed white dwarf
	#with radius R and mass M
	return((R_WD^3/(G*M_WD))^(1/2))
end

# *** #
function Roche_Limit(r,M_WD,M_NS)
	#@assert M_WD < M_NS # In production ready code, consider commenting out assertion
	@assert r > 0 # In production ready code, consider commenting out assertion
	#returns Roche limit of M_WD for companion of mass M_NS 
	#at a separation of distance r 
	q=M_WD/M_NS
    q_2_3rds = q^(2/3)
	return(r * 0.49*q_2_3rds/(0.6*q_2_3rds+log(1+q^(1/3))))
end


# Orbital period
function period(Mtot, a)
	@assert Mtot > 0 # In production ready code, consider commenting out assertion
	@assert a > 0 # In production ready code, consider commenting out assertion
    #return sqrt(4.0 * π^2 * a^3 / (G * Mtot))
    return period_factor * sqrt(a^3 / Mtot)
end


# Decretion rate
function decretion_rate(a, r, R_WD,
                        M_WD, M_NS; A=10.0)
	#@assert R_WD > 0 # In production ready code, consider commenting out assertion
    RL = Roche_Limit(r, M_WD, M_NS)
    if RL > R_WD
        return 0.0
    end
    return -A * M_WD / period(M_WD + M_NS, a) * ((R_WD - RL) / R_WD)^3
end


# Gamma (specific angular momentum of lost matter)
function gamma(M_WD, M_NS; 
               a_ring_frac=1.0, mode::String="Isotropic Re-emission") 
	#@assert M_WD > 0 # In production ready code, consider commenting out assertion
	@assert M_NS > 0 # In production ready code, consider commenting out assertion
	@assert a_ring_frac > 0 # In production ready code, consider commenting out assertion
    first_char = lowercase(first(mode))
    if first_char == 'j'
        return M_NS / M_WD
    elseif first_char == 'i'
        return M_WD / M_NS
    elseif first_char == 'c'
        return (M_WD + M_NS)^2 / (M_WD * M_NS) * sqrt(a_ring_frac)
    else
        error("$mode not defined. Use 'Jeans Mode', 'Isotropic Re-emission', or 'Circumbinary Ring'")
    end
end


# Eddington accretion rate
function eddington_rate(M; mu_e=2.0, eta_acc=0.1)
    @assert M > 0 # In production ready code, consider commenting out assertion
	@assert mu_e > 0 # In production ready code, consider commenting out assertion
	@assert 0 < eta_acc <= 1 # In production ready code, consider commenting out assertion
	#return 1.41e-15 * M * (mu_e / 2.0) * (0.1 / eta_acc)
    return 0.0705e-15 * M * (mu_e / eta_acc) 
end


# Beta (fraction of accreted mass)
function beta(a, r, R_WD, M_WD, M_NS; 
              A=10.0, mu_e=2.0, eta_acc=0.1)
    M_WD_dot = -decretion_rate(a, r, R_WD, M_WD, M_NS; A=A)
    lam_edd = eddington_rate(M_NS; mu_e=mu_e, eta_acc=eta_acc)
    if M_WD_dot < lam_edd
        return 1.0
    end
    return lam_edd / M_WD_dot
end


# Angular momentum loss rate
function J_dot(J, a, r, R_WD, M_WD, M_NS; 
               A=10.0, mu_e=2.0, eta_acc=0.1, 
               a_ring_frac=1.0, mode::String="Isotropic")
    g = gamma(M_WD, M_NS; a_ring_frac=a_ring_frac, mode=mode)
    b = beta(a, r, R_WD, M_WD, M_NS; A=A, mu_e=mu_e, eta_acc=eta_acc)
	@assert b <= 1.0 # In production-ready code, consider commenting this out
    M_WD_dot = decretion_rate(a, r, R_WD, M_WD, M_NS; A=A)
    return J * g * (1.0 - b) * M_WD_dot / (M_WD + M_NS)
end

# *** #
# Semi-major axis change rate
function a_dot(a, r, R_WD, M_WD, M_NS; 
               A=10.0, mu_e=2.0, eta_acc=0.1, 
               a_ring_frac=1.0, mode::String="Isotropic") 
    # Promote keyword values to T so arithmetic mixes consistently
    g = gamma(M_WD, M_NS; a_ring_frac=a_ring_frac, mode=mode)
    b = beta(a, r, R_WD, M_WD, M_NS; A=A, mu_e=mu_e, eta_acc=eta_acc)
	@assert b <= 1.0 # In production-ready code, consider commenting this out
    M_WD_dot = decretion_rate(a, r, R_WD, M_WD, M_NS; A=A)
    one_minus_b = 1 - b
    return -2.0 * a * M_WD_dot / M_WD * (one_minus_b * M_WD / M_NS - one_minus_b * (g + 0.5) * M_WD / (M_WD + M_NS))
end


# White dwarf radius change
function dot_R_WD(M_WD, R_WD)
    natural_R = R0_WD(M_WD)
    tau = dynamical_timescale(M_WD, R_WD)
    return (natural_R - R_WD) / tau
end


# Reduced mass
function mu(M_WD, M_NS)
	#@assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    return M_WD * M_NS / (M_WD + M_NS)
end

#=
# Angular velocity derivative
function dot_phi(J::T, r::T, M_WD::T, M_NS::T)::T where T<:Real
	@assert r > 0 # In production-ready code, consider commenting this out
    return J / (mu(M_WD,M_NS) * r)
end


# Unit vector in the direction of x
function n_hat(x::Vector{Float64})::Vector{Float64}
    return x / norm(x)
end


# Clockwise unit vector perpendicular to x
function theta_hat(x::Vector{Float64})::Vector{Float64}
    n::Vector{Float64} = n_hat(x)
	@assert norm(n) ≈ 1.0 # In production-ready code, consider commenting this out
    return [-n[2], n[1]]
end


# Time derivative of reduced mass
function dot_mu(M_WD::Float64, M_NS::Float64, dM_WD::Float64, dM_NS::Float64)::Float64
	@assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    return (dM_WD * M_NS + M_WD * dM_NS) / (M_WD + M_NS) - M_WD * M_NS * (dM_WD + dM_NS) / (M_WD + M_NS)^2
end


# Symmetric trace-free (STF) operator
function SFT(tensor::Matrix{Float64})::Matrix{Float64}
    sym_part::Matrix{Float64} = (tensor + transpose(tensor)) / 2.0
    tr_val::Float64 = tr(tensor)
    n::Int = size(tensor, 1)
    return sym_part - (1.0 / 3.0) * tr_val * I(n)
end

=#
# Angular momentum loss due to radiation reaction
function dJ_rr(I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
    # Corresponds to: -2G/(5c^5) * (dot(I2[0],I3[1]) - dot(I2[1],I3[0]))
    term::Float64 = dot(I2[:, 1], I3[:, 2]) - dot(I2[:, 2], I3[:, 1])
    #return -2.0 * G / (5.0 * c^5) * term
    return dJ_rr_factor * term
end

#=
# Binary hardening due to radiation reaction
function da_rr(a::Float64, J::Float64, I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
	@assert a > 0 # In production-ready code, consider commenting this out
    return 2.0 * a * dJ_rr(I2, I3) / J
end


# Derivative of μ²
function dmu_2(mu::Float64, dmu::Float64)::Float64
	@assert mu > 0 # In production-ready code, consider commenting this out
    return 2.0 * mu * dmu
end


begin
	# First derivative of θ
	function dtheta(J::Float64, r::Float64, mu::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
		@assert mu > 0 # In production-ready code, consider commenting this out
	    return J / (r * mu)
	end
	
	# Second derivative of θ
	function ddtheta(r::Float64, mu::Float64, dtheta::Float64,
	                 dmu::Float64, dr::Float64, dJ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
		@assert mu > 0 # In production-ready code, consider commenting this out
	    num::Float64 = dJ - dr * mu * dtheta - r * dmu * dtheta
	    den::Float64 = r * mu
	    return num / den
	end
end


begin
# Polar coordinate helpers=#


function x(r, θ)
    @assert r > 0 # In production-ready code, consider commenting this out
    return r * cos(θ)
end

function y(r, θ)
    @assert r > 0 # In production-ready code, consider commenting this out
    return r * sin(θ)
end
	
	#=function dx(r::Float64, θ::Float64, dr::Float64, dθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return dr * cos(θ) - y(r, θ) * dθ
	end
	
	function dy(r::Float64, θ::Float64, dr::Float64, dθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return dr * sin(θ) + x(r, θ) * dθ
	end
	
	function ddx(r::Float64, θ::Float64, dr::Float64, dθ::Float64,
	             ddr::Float64, ddθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return ddr * cos(θ) - 2.0 * dr * sin(θ) * dθ -
	           x(r, θ) * dθ^2 - y(r, θ) * ddθ
	end
	
	function ddy(r::Float64, θ::Float64, dr::Float64, dθ::Float64,
	             ddr::Float64, ddθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return ddr * sin(θ) + 2.0 * dr * cos(θ) * dθ -
	           y(r, θ) * dθ^2 + x(r, θ) * ddθ
	end
end


# Quadrupole moment for effective one-body system
function I_qp(M_WD::Float64, M_NS::Float64, x::Float64, y::Float64)::Matrix{Float64}
	@assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    t_mu::Float64 = mu(M_WD, M_NS)
    return t_mu * [x^2  x * y;
                   x * y  y^2]
end
=#