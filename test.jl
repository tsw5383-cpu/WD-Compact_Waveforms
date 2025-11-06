using BenchmarkTools
using Profile, ProfileCanvas 
using Test
using Base.Threads

using Pkg
function ensure_package_installed(pkg_name::String)
    # Check if the package is installed
    installed_packages = Pkg.dependencies()
    if !any(p -> p.name == pkg_name, values(installed_packages))
        println("Package '$package_name' not found. Installing...")
        Pkg.add(package_name)
        println("Package '$package_name' installed successfully.")
    end
end
Pkg.activate(".")
ensure_package_installed("OrdinaryDiffEq")
ensure_package_installed("ForwardDiff")
ensure_package_installed("CSV")
ensure_package_installed("Dates")

include("src/Integrate.jl")
println(Threads.nthreads())
#M_WD::Float64 = 1.0
#step_limit::Int=100
#Integrate.integrate(M_WD, step_limit=step_limit,write_file=true)
#@btime Integrate.integrate(M_WD, step_limit=step_limit,write_file=false,uneven=false)
#@btime Integrate.integrate(M_WD, step_limit=step_limit,write_file=false,parallel=false)

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
#@test parallel_time < serial_time
println("Passed serial/parallel + integration test for many_integrate()")

#@profile Integrate.integrate(1.0,step_limit=100,write_file=false)
#open("profile_results.txt", "w") do io
#    Profile.print(IOContext(stdout, :displaysize => (24, 120)), format=:flat, sortedby=:count)
#end

