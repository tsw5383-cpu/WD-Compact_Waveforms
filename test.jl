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
ensure_package_installed("OrdinaryDiffEq")
ensure_package_installed("ForwardDiff")
ensure_package_installed("CSV")
ensure_package_installed("Dates")

include("src/integrate.jl")
println(Threads.nthreads())
#M_WD::Float64 = 1.0
#step_limit::Int=100
#Integrate.integrate(M_WD, step_limit=step_limit,write_file=true)
#@btime Integrate.integrate(M_WD, step_limit=step_limit,write_file=false)


@profile Integrate.integrate(1.0,step_limit=100,write_file=false)
open("profile_results.txt", "w") do io
    Profile.print(IOContext(stdout, :displaysize => (24, 120)), format=:flat, sortedby=:count)
end

