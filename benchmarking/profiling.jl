### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ cfb0a424-d4a4-11f0-2222-892b0bd19b98
begin
	using Profile
	using Pkg
	using ProfileCanvas
	Pkg.activate(".")
	using Integrate
end

# ╔═╡ c0956607-797f-4060-9a0e-69f9e6203d9b
let
	Integrate.integrate(1.0; step_limit=2)  # make sure its compiled
	ProfileCanvas.@profview Integrate.integrate(1.0; step_limit=20)
end

# ╔═╡ 166b645e-4e83-4119-a471-6a5bea73cba1
begin
	Profile.clear()   # Clear previous profiling data
	Integrate.integrate(1.0; step_limit=2)
	@profile begin
	    # Call the code you want to profile
	    Integrate.integrate(1.0; step_limit=100)
	end
end

# ╔═╡ 31458fa7-944f-48e0-a38a-d63ab98e42d7
Profile.print(IOContext(stdout, :displaysize => (24, 120)), format=:flat, sortedby=:count)

# ╔═╡ Cell order:
# ╠═cfb0a424-d4a4-11f0-2222-892b0bd19b98
# ╠═c0956607-797f-4060-9a0e-69f9e6203d9b
# ╠═166b645e-4e83-4119-a471-6a5bea73cba1
# ╠═31458fa7-944f-48e0-a38a-d63ab98e42d7
