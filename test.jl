using BenchmarkTools
using Profile, ProfileCanvas 
using Test

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

M_WD::Float64 = 1.0
step_limit::Int=100
@time Integrate.integrate(M_WD, step_limit=step_limit,write_file=true)

#using Profile
#@profile integrate(1.0,step_limit=10000)
#Profile.print()
