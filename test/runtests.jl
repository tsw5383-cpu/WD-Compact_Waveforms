using Test

@testset "Running victor.jl" begin
    include("../victor.jl")
end;

# Type tests
@testset "Type tests" begin 
    M_WD::Float64 = 1.0
    M_NS::Float64 = 1.7
    M_tot::Float64 = M_WD + M_NS
    r::Float64 = 2.2
    a::Float64 = 2.3
    R_WD::Float64 = 1.0
    J::Float64 = 3.3
    @test isa(R_WD(M_WD), Float64)
    @test isa(dynamical_timescale(M_WD, R_WD(M_WD)), Float64)
    @test isa(Roche_Limit(r, M_WD, M_NS), Float64)
    @test isa(period(M_tot, a), Float64)
    @test isa(decretion_rate(a, r, R_WD, M_WD, M_NS), Float64)
    @test isa(gamma(M_WD, M_NS, mode=="j"), Float64)
    @test isa(gamma(M_WD, M_NS, mode=="i"), Float64)
    @test isa(gamma(M_WD, M_NS, mode=="c"), Float64)
    @test isa(eddington_rate(M_NS), Float64)
    @test isa(beta(a, r, R_WD, M_WD, M_NS), Float64)
    @test isa(Jdot(J, a, r, R_WD, M_WD, M_NS), Float64)
    @test isa(a_dot(a, r, R_WD, M_WD, M_NS), Float64)
    @test isa(dot_R_WD(M_WD, R_WD), Float64)
    @test isa(mu(M_WD, M_NS), Float64)
    @test isa(dot_phi(J, r, M_WD, M_NS), Float64)
    x::Vector{Float64} = ones(3) *. 1.63
    @test isa(n_hat(x), Float64)
    @test isa(theta_hat(x), Float64)
    dM_WD::Float64 = -0.4
    dM_NS::Float64 = 0.34
    @test isa(dot_mu(M_WD, M_NS, dM_WD, dM_NS), Float64)
    I2::Matrix{Float64} = Float64[1.0 2.0; 3.0 4.0]
    I3::Matrix{Float64} = Float64[1.2 2.3; -3.1 4.3]
    @test isa(dJ_rr(I2, I3), Float64)
    @test isa(da_rr(a, J, I2, I3), Float64)
    mu::Float64 = 2.0
    dmu::Float64 = 0.2
    @test isa(dmu_2(mu, dmu), Float64)
    @test isa(dtheta(J, r, mu), Float64)
    dr::Float64 = -0.03
    dJ::Float64 = -2.0
    @test isa(ddtheta(r, mu, dtheta(J, r, mu), dmu, dr, dJ), Float64)
    theta::Float64 = 1.03
    @test isa(x(r, theta), Float64)
    @test isa(y(r, theta), Float64)
    @test isa(dx(r, theta, dr, dtheta(J, r, mu)), Float64)
    @test isa(dy(r, theta, dr, dtheta(J, r, mu)), Float64)
    ddr::Float64 = -0.01
    ddtheta::Float64 = 0.02
    @test isa(ddx(r, theta, dr, dtheta(J, r, mu), ddr, ddtheta), Float64)
    @test isa(ddy(r, theta, dr, dtheta(J, r, mu), ddr, ddtheta), Float64)
    @test isa(I(M_WD, M_NS, x(r, theta), y(r, theta)), Matrix{Float64})
    tensor::Matrix{Float64} = Float64[1.0 2.2; 3.4 5.4]
    @test isa(SFT(tensor), Matrix{Float64})
end;

# Reasonable parameter (assertion) tests
@testset "Assertion tests" begin
    @test_throws AssertionError R_WD(M_WD=Float64(-1.0))
    @test_throws AssertionError R_WD(M_WD=Float64(1.0), mu_e=Float64(-1.0))
    @test_throws AssertionError dynamical_timescale(M_WD=Float64(-1.0), R_WD=Float64(1.0))
    @test_throws AssertionError dynamical_timescale(M_WD=Float64(1.0), R_WD=Float64(-1.0))
    @test_throws AssertionError Roche_Limit(r=Float64(-1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
    @test_throws AssertionError Roche_Limit(r=Float64(1.0), M_WD=Float64(-1.0), M_NS=Float64(1.0))
    @test_throws AssertionError Roche_Limit(r=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(-1.0))
    @test_throws AssertionError period(M_tot=Float64(-1.0), a=Float64(1.0))
    @test_throws AssertionError period(M_tot=Float64(1.0), a=Float64(-1.0))
    @test_throws AssertionError decretion_rate(a=Float64(-1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
    @test_throws AssertionError decretion_rate(a=Float64(1.0), r=Float64(-1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
    @test_throws AssertionError decretion_rate(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(-1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
    @test_throws AssertionError decretion_rate(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(-1.0), M_NS=Float64(1.0))
    @test_throws AssertionError decretion_rate(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(-1.0))
    @test_throws AssertionError decretion_rate(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), A=Float64(-1.0))
    @test_throws AssertionError gamma(M_WD=Float64(-1.0), M_NS=Float64(1.0), mode="j")
    @test_throws AssertionError gamma(M_WD=Float64(1.0), M_NS=Float64(-1.0), mode="j")
    @test_throws AssertionError gamma(M_WD=Float64(1.0), M_NS=Float64(1.0), a_ring_frac=Float64(-1.0), mode="j")
    @test_throws AssertionError eddington_rate(M_NS=Float64(-1.0))
    @test_throws AssertionError eddington_rate(M_NS=Float64(1.0), mu_e=Float64(-1.0))
    @test_throws AssertionError eddington_rate(M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(-1.0))
    @test_throws AssertionError eddington_rate(M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(2.0))
      @test_throws AssertionError beta(a=Float64(-1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(-1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(-1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(-1.0), M_NS=Float64(1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(-1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(-1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(-1.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(2.0))
      @test_throws AssertionError beta(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), A=Float64(-1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(-1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(-1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(-1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(-1.0), M_NS=Float64(1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(-1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(-1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(-1.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(2.0))
      @test_throws AssertionError Jdot(J=Float64(1.0), a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), A=Float64(-1.0))
      @test_throws AssertionError a_dot(a=Float64(-1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(-1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(-1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(-1.0), M_NS=Float64(1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(-1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(-1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(-1.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), mu_e=Float64(1.0), eta_acc=Float64(2.0))
      @test_throws AssertionError a_dot(a=Float64(1.0), r=Float64(1.0), R_WD=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(1.0), A=Float64(-1.0))
      @test_throws AssertionError dot_R_WD(M_WD=Float64(-1.0), R_WD=Float64(1.0))
      @test_throws AssertionError dot_R_WD(M_WD=Float64(1.0), R_WD=Float64(-1.0))
      @test_throws AssertionError mu(M_WD=Float64(-1.0), M_NS=Float64(1.0))
      @test_throws AssertionError mu(M_WD=Float64(1.0), M_NS=Float64(-1.0))
      @test_throws AssertionError dot_phi(J=Float64(1.0), r=Float64(-1.0), M_WD=Float64(1.0), M_NS=Float64(1.0))
      @test_throws AssertionError dot_phi(J=Float64(1.0), r=Float64(1.0), M_WD=Float64(-1.0), M_NS=Float64(1.0))
      @test_throws AssertionError dot_phi(J=Float64(1.0), r=Float64(1.0), M_WD=Float64(1.0), M_NS=Float64(-1.0))
      @test_throws AssertionError dot_mu(M_WD=Float64(-1.0), M_NS=Float64(1.0), dM_WD=Float64(1.0), dM_NS=Float64(1.0))
      @test_throws AssertionError dot_mu(M_WD=Float64(1.0), M_NS=Float64(-1.0), dM_WD=Float64(1.0), dM_NS=Float64(1.0))
      I2::Matrix{Float64} = Float64[1.0 2.0; 3.0 4.0]
      I3::Matrix{Float64} = Float64[1.2 2.3; -3.1 4.3]
      @test_throws AssertionError da_rr(a=Float64(-1.0), J=Float64(1.0), I2=I2, I3=I3)
      @test_throws AssertionError dmu_2(mu=Float64(-1.0), dmu=Float64(1.0))
      @test_throws AssertionError dtheta(J=Float64(1.0), r=Float64(-1.0), mu=Float64(1.0))
      @test_throws AssertionError dtheta(J=Float64(1.0), r=Float64(1.0), mu=Float64(-1.0))
      @test_throws AssertionError ddtheta(r=Float64(-1.0), mu=Float64(1.0), dtheta=Float64(1.0), dmu=Float64(1.0), dr=Float64(1.0), dJ=Float64(1.0))
      @test_throws AssertionError ddtheta(r=Float64(1.0), mu=Float64(-1.0), dtheta=Float64(1.0), dmu=Float64(1.0), dr=Float64(1.0), dJ=Float64(1.0))
      @test_throws AssertionError x(r=Float64(-1.0), theta=Float64(1.0))
      @test_throws AssertionError y(r=Float64(-1.0), theta=Float64(1.0))
      @test_throws AssertionError dx(r=Float64(-1.0), theta=Float64(1.0), dr=Float64(1.0), dtheta=Float64(1.0))
      @test_throws AssertionError dy(r=Float64(-1.0), theta=Float64(1.0), dr=Float64(1.0), dtheta=Float64(1.0))
      @test_throws AssertionError ddx(r=Float64(-1.0), theta=Float64(1.0), dr=Float64(1.0), dtheta=Float64(1.0), ddr=Float64(1.0), ddtheta=Float64(1.0))
      @test_throws AssertionError ddy(r=Float64(-1.0), theta=Float64(1.0), dr=Float64(1.0), dtheta=Float64(1.0), ddr=Float64(1.0), ddtheta=Float64(1.0))
      @test_throws AssertionError I(M_WD=Float64(-1.0), M_NS=Float64(1.0), x=Float64(1.0), y=Float64(1.0))
      @test_throws AssertionError I(M_WD=Float64(1.0), M_NS=Float64(-1.0), x=Float64(1.0), y=Float64(1.0))
end;

# Reasonable input -> output mapping tests
@testset "Expected output tests" begin
   @test R0_WD(M_WD=1e-10) > 0
   @test dynamical_timescale(M_WD=1e-10, R_WD=1e-10) > 0
   @test Roche_Limit(r=1e-10, M_WD=1e-10, M_NS=1e-10) > 0
   @test period(M_tot=1e-10, a=1e-10) > 0
   @test decretion_rate(a=1e-10, r=1e-10, R_WD=1e-20, M_WD=1e-10, M_NS=1e-10) == 0.0 # Here, Roche limit is bigger than R_WD
   @test gamma(M_WD=1e-10, M_NS=1e-10, mode="j") > 0
   @test gamma(M_WD=1e-10, M_NS=1e-10, mode="i") > 0
   @test gamma(M_WD=1e-10, M_NS=1e-10, mode="c") > 0
   @test eddington_rate(M_NS=1e-10) > 0
   @test beta(a=1e-10, r=1e-10, R_WD=1e-20, M_WD=1e-10, M_NS=1e-10) == 1.0 # Here, Roche limit is bigger than R_WD, so M_dot = 0.0, so beta = 1.0
   @test 0.0 < beta(a=1e-10, r=1e-10, R_WD=1e-10, M_WD=1e-10, M_NS=1e-10) < 1.0 # Here, Roche limit is not than R_WD, so M_dot > 0
   @test mu(M_WD=1e-10, M_NS=1e-10) > 0
   x::Vector{Float64} = ones(3) .* 1e-10
   @test norm(n_hat(x)) ~ 1.0, rtol=1e-5
   @test norm(theta_hat(x)) ~ 1.0, rtol=1e-5
end;

# julia-actions/julia-runtest
touch(joinpath(ENV["HOME"], "julia-runtest"))
