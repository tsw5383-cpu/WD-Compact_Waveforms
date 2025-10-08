import Pkg
Pkg.add("BenchmarkTools")
Pkg.add("ProfileCanvas")
Pkg.add("PlutoUI")
Pkg.add("PlutoTest")
Pkg.add("PlutoTeachingTools")
Pkg.add("PlutoLinks")
Pkg.add("Test")
Pkg.add("LinearAlgebra")

using Test
using LinearAlgebra

@testset "Running tests" begin
    
    include("../victor.jl")
	    	
	# Type tests
	@testset "Type tests" begin 
	    M_WD::Float64 = 1.0
	    M_NS::Float64 = 1.7
	    M_tot::Float64 = M_WD + M_NS
	    r::Float64 = 2.2
	    a::Float64 = 2.3
	    J::Float64 = 3.3
		R_WD_test = R_WD(M_WD)
	    @test isa(R_WD_test, Float64)
	    @test isa(dynamical_timescale(M_WD, R_WD_test), Float64)
	    @test isa(Roche_Limit(r, M_WD, M_NS), Float64)
	    @test isa(period(M_tot, a), Float64)
	    @test isa(decretion_rate(a, r, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(gamma(M_WD, M_NS, mode="j"), Float64)
	    @test isa(gamma(M_WD, M_NS, mode="i"), Float64)
	    @test isa(gamma(M_WD, M_NS, mode="c"), Float64)
	    @test isa(eddington_rate(M_NS), Float64)
	    @test isa(beta(a, r, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(J_dot(J, a, r, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(a_dot(a, r, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(dot_R_WD(M_WD, R_WD_test), Float64)
		mu_test = mu(M_WD, M_NS)
	    @test isa(mu_test, Float64)
	    @test isa(dot_phi(J, r, M_WD, M_NS), Float64)
	    x_vec::Vector{Float64} = ones(2) .* 1.32
	    @test isa(n_hat(x_vec), Vector{Float64})
	    @test isa(theta_hat(x_vec), Vector{Float64})
	    dM_WD::Float64 = -0.4
	    dM_NS::Float64 = 0.34
	    @test isa(dot_mu(M_WD, M_NS, dM_WD, dM_NS), Float64)
	    I2::Matrix{Float64} = Float64[1.0 2.0; 3.0 4.0]
	    I3::Matrix{Float64} = Float64[1.2 2.3; -3.1 4.3]
	    @test isa(dJ_rr(I2, I3), Float64)
	    @test isa(da_rr(a, J, I2, I3), Float64)
	    dmu::Float64 = 0.2
	    @test isa(dmu_2(mu_test, dmu), Float64)
		dtheta_test = dtheta(J, r, mu_test)
	    @test isa(dtheta_test, Float64)
	    dr::Float64 = -0.03
	    dJ::Float64 = -2.0
		ddtheta_test = ddtheta(r, mu_test, dtheta_test, dmu, dr, dJ)
	    @test isa(ddtheta_test, Float64)
	    theta::Float64 = 1.03
		x_test = x(r, theta)
		y_test = y(r, theta)
	    @test isa(x_test, Float64)
	    @test isa(y_test, Float64)
	    @test isa(dx(r, theta, dr, dtheta_test), Float64)
	    @test isa(dy(r, theta, dr, dtheta_test), Float64)
	    ddr::Float64 = -0.01
	    @test isa(ddx(r, theta, dr, dtheta_test, ddr, ddtheta_test), Float64)
	    @test isa(ddy(r, theta, dr, dtheta_test, ddr, ddtheta_test), Float64)
	    #@test isa(I_qp(M_WD, M_NS, x_test, y_test), Matrix{Float64})
	    tensor::Matrix{Float64} = [1.0 2.2; 3.4 5.4]
	    @test isa(SFT(tensor), Matrix{Float64})
	end;
	
	# Reasonable input -> output mapping tests
	@testset "Expected output tests" begin
	   @test R_WD(1e-10) > 0
	   @test dynamical_timescale(1e-10, 1e-10) > 0
	   @test Roche_Limit(1e-10, 1e-10, 1e-10) > 0
	   @test period(1e-10, 1e-10) > 0
	   @test decretion_rate(1e-10, 1e-10, 1e-20, 1e-10, 1e-10) == 0.0 # Here, Roche limit is bigger than R_WD
	   @test gamma(1e-10, 1e-10, mode="j") > 0
	   @test gamma(1e-10, 1e-10, mode="i") > 0
	   @test gamma(1e-10, 1e-10, mode="c") > 0
	   @test eddington_rate(1e-10) > 0
	   @test beta(1e-10, 1e-10, 1e-20, 1e-10, 1e-10) == 1.0 # Here, Roche limit is bigger than R_WD, so M_dot = 0.0, so beta = 1.0
	   @test 0.0 < beta(1e-10, 1e-10, 1e-10, 1e-10, 1e-10) < 1.0 # Here, Roche limit is not than R_WD, so M_dot > 0
	   @test mu(1e-10, 1e-10) > 0
		x_vec::Vector{Float64} = ones(2) .* 1e-10
	   @test norm(n_hat(x_vec)) ≈ 1.0 rtol=1e-5
	   @test norm(theta_hat(x_vec)) ≈ 1.0 rtol=1e-5
	end;

end;

# julia-actions/julia-runtest
touch(joinpath(ENV["HOME"], "julia-runtest"))
