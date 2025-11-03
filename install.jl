using Pkg

"""Checks if the package is installed. If not, installs the package."""
function ensure_package_installed(pkg_name::String)
    installed_packages = Pkg.dependencies()
    if !any(p -> p.name == pkg_name, values(installed_packages))
        println("Package '$package_name' not found. Installing...")
        Pkg.add(package_name)
        println("Package '$package_name' installed successfully.")
    end
end

ensure_package_installed("Memoization")
ensure_package_installed("LinearAlgebra")
ensure_package_installed("OrdinaryDiffEq")
ensure_package_installed("ForwardDiff")
ensure_package_installed("CSV")
ensure_package_installed("Dates")

Pkg.precompile()