using Pkg
Pkg.add("BenchmarkTools")
Pkg.add("PlutoTest")
using Test
using BenchmarkTools

@testset "Running tests" begin
    #using Pkg
    using BenchmarkTools
    using PlutoTest
    using LinearAlgebra
    include("../src/Integrate.jl")
    include("../src/math.jl")
	    	
	# Type tests
	@testset "Type tests" begin 
        println("Beginning type tests")
	    M_WD::Float64 = 1.0
	    M_NS::Float64 = 1.7
	    M_tot::Float64 = M_WD + M_NS
	    a::Float64 = 2.3
	    J::Float64 = 3.3
		R_WD_test = Math.R0_WD(M_WD)
	    @test isa(R_WD_test, Float64)
	    @test isa(Math.dynamical_timescale(M_WD, R_WD_test), Float64)
	    @test isa(Math.Roche_Limit(a, M_WD, M_NS), Float64)
	    @test isa(Math.period(a, M_tot), Float64)
	    @test isa(Math.decretion_rate(a, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(Math.gamma(M_WD, M_NS, mode="j"), Float64)
	    @test isa(Math.gamma(M_WD, M_NS, mode="i"), Float64)
	    @test isa(Math.gamma(M_WD, M_NS, mode="c"), Float64)
	    @test isa(Math.eddington_rate(M_NS), Float64)
	    @test isa(Math.beta(a, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(Math.J_dot(J, a, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(Math.a_dot(a, R_WD_test, M_WD, M_NS), Float64)
	    @test isa(Math.dot_R_WD(M_WD, R_WD_test), Float64)
		mu_test = Math.mu(M_WD, M_NS)
	    @test isa(mu_test, Float64)
	    I2::Matrix{Float64} = Float64[1.0 2.0; 3.0 4.0]
	    I3::Matrix{Float64} = Float64[1.2 2.3; -3.1 4.3]
	    @test isa(Math.dJ_rr(I2, I3), Float64)
        r::Float64 = 90.0
	    theta::Float64 = 1.03
		x_test = Math.x(r, theta)
		y_test = Math.y(r, theta)
	    @test isa(x_test, Float64)
	    @test isa(y_test, Float64)
        println("Passed all type tests")
	end;
	
	# Reasonable input -> output mapping tests
	@testset "Expected output tests" begin
        println("Beginning expected output tests")
	   @test Math.R0_WD(1e-10) > 0
	   @test Math.dynamical_timescale(1e-10, 1e-10) > 0
	   @test Math.Roche_Limit(1e-10, 1e-10, 1e-10) > 0
	   @test Math.period(1e-10, 1e-10) > 0
	   @test Math.decretion_rate(1e-10, 1e-20, 1e-10, 1e-10) == 0.0 # Here, Roche limit is bigger than R_WD
	   @test Math.gamma(1e-10, 1e-10, mode="j") > 0
	   @test Math.gamma(1e-10, 1e-10, mode="i") > 0
	   @test Math.gamma(1e-10, 1e-10, mode="c") > 0
	   @test Math.eddington_rate(1e-10) > 0
	   @test Math.beta(1e-10, 1e-20, 1e-10, 1e-10) == 1.0 # Here, Roche limit is bigger than R_WD, so M_dot = 0.0, so beta = 1.0
	   @test 0.0 < Math.beta(1e-10, 1e-10, 1e-10, 1e-10) < 1.0 # Here, Roche limit is not than R_WD, so M_dot > 0
	   @test Math.mu(1e-10, 1e-10) > 0
        println("Passed all expected output tests")
	end;

    # Comparing serial and parallel versions of code, which also serves as an integration test because testset would fail if integrate() didn't work
    @testset "Serial/parallel test for integrate()" begin
        using BenchmarkTools
        nthreads::Int = Threads.nthreads()
        if nthreads < 2
            throw(error("ERROR: Currently you are running this test with $(Threads.nthreads()) threads. This test requires at least 2 threads. Request more cores before proceeding."))
        end
        println("Beginning serial/parallel + integration test for integrate()")
        println("The number of threads being used in parallel version is $(Threads.nthreads()).")
        M_WD::Float64 = 1.0
        step_limit::Int=100
        parallel_time = @belapsed Integrate.integrate($M_WD, step_limit=$step_limit,write_file=false,uneven=false)
        serial_time = @belapsed Integrate.integrate($M_WD, step_limit=$step_limit,write_file=false,parallel=false)
        println("M_WD: $(M_WD), step_limit: $(step_limit)")
        println("The time for integrate() in the parallel version is $(parallel_time) sec.")
        println("The time for integrate() in the serial version is $(serial_time) sec.")
        @test parallel_time < serial_time
        @test parallel_time < 0.4
        @test serial_time < 0.4 # I got speeds of ~200 ms for serial
        println("Passed serial/parallel + integration test for integrate()")
    end

    @testset "Serial/parallel test for many_integrate()" begin
        sim_conditions = (M_WD_lower=0.5, M_WD_upper=0.8, M_WD_step=0.1, 
                          M_NS_lower=1.2, M_NS_upper=1.4, M_NS_step=0.1)
        step_limit::Int=100
        println("Beginning serial/parallel + integration test for many_integrate()")
        println("The number of threads being used in the parallel version is $(Threads.nthreads()).")
        println("The range of M_WDs: (lower=$(sim_conditions.M_WD_lower),upper=$(sim_conditions.M_WD_upper),step=$(sim_conditions.M_WD_step)). The range of M_NSs: (lower=$(sim_conditions.M_NS_lower),upper=$(sim_conditions.M_NS_upper),step=$(sim_conditions.M_NS_step)). step_limit: $(step_limit)")
        parallel_time_many = @belapsed Integrate.many_integrate($sim_conditions, step_limit=$step_limit,write_file=false,uneven=false)
        serial_time_many = @belapsed Integrate.many_integrate($sim_conditions, step_limit=$step_limit,write_file=false,parallel=false)
        
        println("The time for many_integrate() in the parallel version is $(parallel_time_many) sec.")
        println("The time for many_integrate() in the serial version is $(serial_time_many) sec.")
        @test parallel_time < serial_time
        println("Passed serial/parallel + integration test for many_integrate()")
    end
end;

# julia-actions/julia-runtest
touch(joinpath(ENV["HOME"], "julia-runtest"))
