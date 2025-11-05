using Pkg

"""Checks if the package is installed. If not, installs the package."""
function ensure_package_installed(pkg_name::String)
    installed_packages = Pkg.dependencies()
    if !any(p -> p.name == pkg_name, values(installed_packages))
        println("Package '$pkg_name' not found. Installing...")
        Pkg.add(pkg_name)
        println("Package '$pkg_name' installed successfully.")
    end
end

ensure_package_installed("Memoization")
ensure_package_installed("LinearAlgebra")
ensure_package_installed("OrdinaryDiffEq")
ensure_package_installed("ForwardDiff")
ensure_package_installed("CSV")
ensure_package_installed("Dates")
ensure_package_installed("NFFT")
ensure_package_installed("FFTW")

Pkg.precompile()