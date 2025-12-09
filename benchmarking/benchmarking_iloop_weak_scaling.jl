using Pkg
using BenchmarkTools
Pkg.activate(".")
include("src/Integrate.jl")

step_limit = parse(Int, ARGS[1])
println("Benchmarking performance vs problem size with a fixed number of threads.")
println("We are using $(Threads.nthreads()) threads on one node.")
println("The step limit in this specific run is $(step_limit).")
println("")

M_WD::Float64 = 1.0
println("Running shared memory even parallelization integrate() with default optional arguments and M_WD $(M_WD).")
shared_even_parallel_time = @belapsed Integrate.integrate(M_WD, step_limit=step_limit,write_file=false,uneven=false)
println("The walltime is $(shared_even_parallel_time) seconds.")
println("")

println("Running shared memory uneven parallelization integrate() with default optional arguments and M_WD $(M_WD).")
shared_uneven_parallel_time = @belapsed Integrate.integrate(M_WD, step_limit=step_limit,write_file=false,uneven=true)
println("The walltime is $(shared_uneven_parallel_time) seconds.")
println("")

println("Running serial integrate() with default optional arguments and M_WD $(M_WD).")
serial_time = @belapsed Integrate.integrate(M_WD, step_limit=step_limit,write_file=false,parallel=false)
println("The walltime is $(serial_time) seconds.")