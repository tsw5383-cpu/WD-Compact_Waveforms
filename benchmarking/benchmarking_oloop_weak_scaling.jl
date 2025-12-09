# This script benchmarks strong scaling for the outer loop (many_integrate()) with shared memory parallelization and serial.
using BenchmarkTools
using Distributed
@everywhere include("../src/Integrate.jl")

println("Benchmarking many_integrate() performance vs problem size with fixed number of nodes and cores per node.")
println("We are using $(Threads.nthreads()) workers on this node.")
println("")

sim_conditions = (M_WDs= 0.5:0.1:0.8, M_NSs= 1.2:0.1:1.4)

step_limit::Int = parse(Int, ARGS[1])
println("Inputted M_WDs: $(sim_conditions.M_WDs). Inputted M_NSs: $(sim_conditions.M_NSs). Inputted step_limit: $(step_limit).")
println("")

println("Running distributed memory parallelization many_integrate() with default optional arguments.")
distributed_parallel_time_many = @belapsed Integrate.many_integrate($sim_conditions, step_limit=$step_limit,write_file=false,uneven=false,distributed=true)
println("The time for many_integrate() in the distributed memory parallel version is $(distributed_parallel_time_many) sec.")

println("Running serial many_integrate() with default optional arguments.")
serial_time_many = @belapsed Integrate.many_integrate($sim_conditions, step_limit=$step_limit,write_file=false,parallel=false)
println("The time for many_integrate() in the serial version is $(serial_time_many) sec.")