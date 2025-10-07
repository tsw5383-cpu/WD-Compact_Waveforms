### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 0cac30f1-9101-4ba3-accc-52621bc1d16f
begin
	using LinearAlgebra
	using BenchmarkTools
	using Profile, ProfileCanvas 
	using PlutoUI, PlutoTest, PlutoTeachingTools, PlutoLinks
end;

# ╔═╡ 5a137355-8c63-492d-961c-b0fcb4897eae
md"""
# Victor Liu
"""

# ╔═╡ 47811823-7c56-40b5-805f-8f76b0efab59
md"
## Rubric
"

# ╔═╡ adeb9a16-52c5-4e96-81fa-19fca745ad2b
md"Code performs proposed tasks (1 point)

Comprehensive set of unit tests, at least one integration or regression test (1 point)

Code passes tests (1 point)

Student code uses a version control system effectively (1 point)

Repository includes many regular, small commits (1 point)

Documentation for functions' purpose and design (1 point)

Comprehensive set of assertions (1 point)

Variable/function names consistent, distinctive & meaningful (1 point)

Useful & consistent code formatting & style (1 point)

Code is modular, rather than having chunks of same code copied and pasted (1 point)"

# ╔═╡ 0ee4cf56-090d-4b9b-a266-a592952e2c1e
md"Mass transfer is hard to come up with an integration test for because well there is no known answer! But we can make a best guess for mass transfer and put that into an integration unit test.

For radiation decay and orbits, we know the answer, so these integration tests we can test."

# ╔═╡ 14bfa5a2-aa27-4c05-a1da-d4be8037139d
md"Don't need to optimize code right now, just need to meet rubric. I might try and do some basic optimization and then profile/benchmark later."

# ╔═╡ a615463a-fcaf-4823-8af4-8a04c79428ab
md"
## Importing Packages and Setting Defaults
"

# ╔═╡ 0bad3ea9-d820-485e-8e84-387bf12a3f0e
md"
## Constants
"

# ╔═╡ b313fdb0-3a5a-4999-a235-93cf20b3222a
# Constants
begin
	const G::Float64 = 1.327e11   # km^3 / M⊙ / s^2
	const c::Float64 = 2.998e5    # km/s
end

# ╔═╡ 0f987ec5-7cd0-4370-b72a-2ba22e7f4c09
md"
## Physics
"

# ╔═╡ 3da8b6b5-f5a0-49f3-a03d-36c87584e583
md"
note: if we want to pass in an array of values to these functions (e.g. calculate R_WD for multiple M and mu at once), will need to refactor these functions to allow for vector operations.

TO-DO: add assertions (e.g. mass of white dwarf, radius of white dwarf, q)
"

# ╔═╡ b7a2300a-bc8f-4387-a324-7f3de70d81fc
function R_WD(M_WD::Float64;mu_e=2::Float64)::Float64
	@assert M_WD > 0 # In production ready code, consider commenting out assertion
	@assert mu_e > 0 # In production ready code, consider commenting out assertion
	#white dwarf radius in km as ftn of mass in M_odot and mean molecular weight/electron=2
	return(10^4 * (M_WD/0.7)^(-1/3) * (1-(M_WD/1.4)^(4/3))^(1/2) * (mu_e/2)^(-5/3))
end

# ╔═╡ b7880c1c-dbd2-4bf3-b916-6165c27c7d81
function dynamical_timescale(M_WD::Float64,R_WD::Float64)::Float64
	@assert M_WD > 0 # In production ready code, consider commenting out assertion
	@assert R_WD > 0 # In production ready code, consider commenting out assertion
	#returns dynamical timescale for perturbed white dwarf
	#with radius R and mass M
	return((R_WD^3/(G*M_WD))^(1/2))
end

# ╔═╡ d52feb78-6ddf-4b75-995e-4b2a1e452571
function Roche_Limit(r::Float64,M_WD::Float64,M_WD::Float64)
	#@assert M1 < M2 # In production ready code, consider commenting out assertion
	@assert r > 0 # In production ready code, consider commenting out assertion
	#returns Roche limit of M1 for companion of mass M2 
	#at a separation of distance r 
	q=M_WD/M_WD
	return(r * 0.49*q^(2/3)/(0.6*q^(2/3)+log(1+q^(1/3))))
end

# ╔═╡ 8cb6660f-7181-478f-8c2e-e67c908ad8b7
# Orbital period
function period(Mtot::Float64, a::Float64)::Float64
	@assert Mtot > 0 # In production ready code, consider commenting out assertion
	@assert a > 0 # In production ready code, consider commenting out assertion
    return sqrt(4.0 * π^2 * a^3 / (G * Mtot))
end

# ╔═╡ b5b5c7a1-3108-4a39-929b-f7af515b55c5
# Decretion rate
function decretion_rate(a::Float64, r::Float64, R1::Float64,
                        M1::Float64, M2::Float64; A::Float64=10.0)::Float64
	@assert R1 > 0 # In production ready code, consider commenting out assertion
    RL::Float64 = Roche_Limit(r, M1, M2)
    if RL > R1
        return 0.0
    end
    return -A * M1 / period(M1 + M2, a) * ((R1 - RL) / RL)^3
end

# ╔═╡ 11fa4c68-4fd8-496e-8603-8181dd3940a4
# Gamma (specific angular momentum of lost matter)
function gamma(M1::Float64, M2::Float64; 
               a_ring_frac::Float64=1.0, mode::String="Isotropic Re-emission")::Float64
	@assert M1 > 0 # In production ready code, consider commenting out assertion
	@assert M2 > 0 # In production ready code, consider commenting out assertion
	@assert a_ring_frac > 0 # In production ready code, consider commenting out assertion
    first_char = lowercase(first(mode))
    if first_char == 'j'
        return M2 / M1
    elseif first_char == 'i'
        return M1 / M2
    elseif first_char == 'c'
        return (M1 + M2)^2 / (M1 * M2) * sqrt(a_ring_frac)
    else
        error("$mode not defined. Use 'Jeans Mode', 'Isotropic Re-emission', or 'Circumbinary Ring'")
    end
end

# ╔═╡ 283caef4-4cac-44fd-8fe9-17507f14baf9
# Eddington accretion rate
function eddington_rate(M::Float64; mu_e::Float64=2.0, eta_acc::Float64=0.1)::Float64
    @assert M > 0 # In production ready code, consider commenting out assertion
	@assert mu_e > 0 # In production ready code, consider commenting out assertion
	@assert 0 < eta_acc <= 1 # In production ready code, consider commenting out assertion
	return 1.41e-15 * M * (mu_e / 2.0) * (0.1 / eta_acc)
end

# ╔═╡ 524b7699-c0f2-4904-993c-7943de80e821
# Beta (fraction of accreted mass)
function beta(a::Float64, r::Float64, R1::Float64, M1::Float64, M2::Float64; 
              A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1)::Float64
    M_dot::Float64 = -decretion_rate(a, r, R1, M1, M2; A=A)
    lam_edd::Float64 = eddington_rate(M2; mu_e=mu_e, eta_acc=eta_acc)
    if M_dot < lam_edd
        return 1.0
    end
    return lam_edd / M_dot
end

# ╔═╡ 6c2ec7e6-1868-4424-ac5a-0d3ee97e0979
# Angular momentum loss rate
function J_dot(J::Float64, a::Float64, r::Float64, R1::Float64, M1::Float64, M2::Float64; 
               A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
               a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64
    g::Float64 = gamma(M1, M2; a_ring_frac=a_ring_frac, mode=mode)
    b::Float64 = beta(a, r, R1, M1, M2; A=A, mu_e=mu_e, eta_acc=eta_acc)
	@assert b <= 1.0 # In production-ready code, consider commenting this out
    m1_dot::Float64 = decretion_rate(a, r, R1, M1, M2; A=A)
    return J * g * (1.0 - b) * m1_dot / (M1 + M2)
end

# ╔═╡ 760b9821-5596-46e5-b263-8101d0b2ea42
# Semi-major axis change rate
function a_dot(a::Float64, r::Float64, R1::Float64, M1::Float64, M2::Float64; 
               A::Float64=10.0, mu_e::Float64=2.0, eta_acc::Float64=0.1, 
               a_ring_frac::Float64=1.0, mode::String="Isotropic")::Float64
    g::Float64 = gamma(M1, M2; a_ring_frac=a_ring_frac, mode=mode)
    b::Float64 = beta(a, r, R1, M1, M2; A=A, mu_e=mu_e, eta_acc=eta_acc)
	@assert b <= 1.0 # In production-ready code, consider commenting this out
    m1_dot::Float64 = decretion_rate(a, r, R1, M1, M2; A=A)
    return -2.0 * a * m1_dot / M1 * (1.0 - b * M1 / M2 - (1.0 - b) * (g + 0.5) * M1 / (M1 + M2))
end

# ╔═╡ 65c9b6ea-de2f-4051-9305-6219cc814bd6
# White dwarf radius change
function dot_R_WD(M1::Float64, R1::Float64)::Float64
    natural_R::Float64 = R_WD(M1)
    tau::Float64 = dynamical_timescale(M1, R1)
    return (natural_R - R1) / tau
end

# ╔═╡ e58e91a2-7f10-460d-9abe-155cd98e36b3
# Reduced mass
function mu(M1::Float64, M2::Float64)::Float64
	@assert M1 > 0 # In production-ready code, consider commenting this out
	@assert M2 > 0 # In production-ready code, consider commenting this out
    return M1 * M2 / (M1 + M2)
end

# ╔═╡ d589347d-e0f4-40f7-b79b-f10821e19534
# Angular velocity derivative
function dot_phi(J::Float64, r::Float64, M1::Float64, M2::Float64)::Float64
	@assert r > 0 # In production-ready code, consider commenting this out
    mu::Float64 = mu(M1,M2)
    return J / (mu * r)
end

# ╔═╡ b5f7928e-0299-426c-86c5-bf620bd3da03
# Unit vector in the direction of x
function n_hat(x::Vector{Float64})::Vector{Float64}
    return x / norm(x)
end

# ╔═╡ 4347c02f-f33f-4a9c-a292-15689edafce1
# Clockwise unit vector perpendicular to x
function theta_hat(x::Vector{Float64})::Vector{Float64}
    n::Vector{Float64} = n_hat(x)
	@assert norm(n) ~ 1.0 # In production-ready code, consider commenting this out
    return [-n[2], n[1]]
end

# ╔═╡ 51133208-6ebf-41c4-bf52-7f5d8e294d0c
# Time derivative of reduced mass
function dot_mu(M1::Float64, M2::Float64, dM1::Float64, dM2::Float64)::Float64
	@assert M1 > 0 # In production-ready code, consider commenting this out
	@assert M2 > 0 # In production-ready code, consider commenting this out
    return (dM1 * M2 + M1 * dM2) / (M1 + M2) - M1 * M2 * (dM1 + dM2) / (M1 + M2)^2
end

# ╔═╡ c1fed078-6b90-4258-a4be-87655ec3907c
# Angular momentum loss due to radiation reaction
function dJ_rr(I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
    # Corresponds to: -2G/(5c^5) * (dot(I2[0],I3[1]) - dot(I2[1],I3[0]))
    term::Float64 = dot(I2[:, 1], I3[:, 2]) - dot(I2[:, 2], I3[:, 1])
    return -2.0 * G / (5.0 * c^5) * term
end

# ╔═╡ a26d7950-3b81-4a09-84ab-5e2c8fb4f00c
# Binary hardening due to radiation reaction
function da_rr(a::Float64, J::Float64, I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
	@assert a > 0 # In production-ready code, consider commenting this out
    return 2.0 * a * dJ_rr(I2, I3) / J
end

# ╔═╡ 0cbb83e3-7ec0-4498-9777-d6a3f56a422e
# Derivative of μ²
function dmu_2(mu::Float64, dmu::Float64)::Float64
	@assert mu > 0 # In production-ready code, consider commenting this out
    return 2.0 * mu * dmu
end

# ╔═╡ f49155de-8db9-45b1-84f2-a85ed55a0dc5
begin
	# First derivative of θ
	function dtheta(J::Float64, r::Float64, mu::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
		@assert mu > 0 # In production-ready code, consider commenting this out
	    return J / (r * mu)
	end
	
	# Second derivative of θ
	function ddtheta(r::Float64, mu::Float64, dtheta::Float64,
	                 dmu::Float64, dr::Float64, dJ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
		@assert mu > 0 # In production-ready code, consider commenting this out
	    num::Float64 = dJ - dr * mu * dtheta - r * dmu * dtheta
	    den::Float64 = r * mu
	    return num / den
	end
end

# ╔═╡ 8bd09c6b-fedb-4256-aa10-ef6f76e9397f
begin
	# Polar coordinate helpers
	function x(r::Float64, θ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return r * cos(θ)
	end
	
	function y(r::Float64, θ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return r * sin(θ)
	end
	
	function dx(r::Float64, θ::Float64, dr::Float64, dθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return dr * cos(θ) - y(r, θ) * dθ
	end
	
	function dy(r::Float64, θ::Float64, dr::Float64, dθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return dr * sin(θ) + x(r, θ) * dθ
	end
	
	function ddx(r::Float64, θ::Float64, dr::Float64, dθ::Float64,
	             ddr::Float64, ddθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return ddr * cos(θ) - 2.0 * dr * sin(θ) * dθ -
	           x(r, θ) * dθ^2 - y(r, θ) * ddθ
	end
	
	function ddy(r::Float64, θ::Float64, dr::Float64, dθ::Float64,
	             ddr::Float64, ddθ::Float64)::Float64
		@assert r > 0 # In production-ready code, consider commenting this out
	    return ddr * sin(θ) + 2.0 * dr * cos(θ) * dθ -
	           y(r, θ) * dθ^2 + x(r, θ) * ddθ
	end
end

# ╔═╡ 46e2e3bd-0519-4d0c-a748-fcff17118b6e
# Quadrupole moment for effective one-body system
function I(M_WD::Float64, M_NS::Float64, x::Float64, y::Float64)::Matrix{Float64}
	@assert M_WD > 0 # In production-ready code, consider commenting this out
	@assert M_NS > 0 # In production-ready code, consider commenting this out
    t_mu::Float64 = mu(M_WD, M_NS)
    return t_mu * [x^2  x * y;
                   x * y  y^2]
end

# ╔═╡ a22422ea-683b-4d38-ac2b-d9802c0c2c80
# Symmetric trace-free (STF) operator
function SFT(tensor::Matrix{Float64})::Matrix{Float64}
    sym_part::Matrix{Float64} = (tensor + transpose(tensor)) / 2.0
    tr_val::Float64 = tr(tensor)
    n::Int = size(tensor, 1)
    return sym_part - (1.0 / 3.0) * tr_val * I(n)
end

# ╔═╡ d293006e-7bd4-46dc-82ec-dcabb7ccc5ed
md"
## Tests (may be moved to a separate Julia file later)
"

# ╔═╡ 0d5207be-9594-4d28-8fa5-acd0cd397acf
# ╠═╡ disabled = true
#=╠═╡
begin
	thetas = [1.0, 1.0, 1.4, 0.01, 1.0, 0.0]
    dthetas = [0.1, 1e-5, 1e-6, 0.0, 1e-3, 0.05]
    ddIs = [rand(2,2), rand(2,2)]
    drs = [0.0]
    dM_WDs = [0.0]
    dt = 0.01

    dJ, da, ddI = fixed_RR_forces(thetas, dthetas, ddIs, drs, dM_WDs, dt)

    @test isa(dJ, Float64)
    @test isa(da, Float64)
    @test isa(ddI, Matrix{Float64})
    @test size(ddI) == (2,2)
    @test isfinite(dJ)
    @test isfinite(da)
    @test all(isfinite, ddI)
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoLinks = "0ff47ea0-7a50-410d-8455-4348d5de0420"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
ProfileCanvas = "efd6af41-a80b-495e-886c-e51b0c7d77a3"

[compat]
BenchmarkTools = "~1.6.0"
PlutoLinks = "~0.1.6"
PlutoTeachingTools = "~0.4.5"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.71"
ProfileCanvas = "~0.1.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "1c26fab669920113f1f3b547670887b14ba9c85a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "e38fbc49a620f5d0b660d7f543db1009fe0f8336"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.6.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "5ac098a7c8660e217ffac31dc2af0964a8c3182a"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "d8337622fe53c05d16f031df24daf0270e53bc64"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.5"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "52e1296ebbde0db845b356abbbe67fb82a0a116c"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.9"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "73b98709ad811a6f81d84e105f4f695c229385ba"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.3"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "85778cdf2bed372008e6646c64340460764a5b85"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.5"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.ProfileCanvas]]
deps = ["Base64", "JSON", "Pkg", "Profile", "REPL"]
git-tree-sha1 = "e42571ce9a614c2fbebcaa8aab23bbf8865c624e"
uuid = "efd6af41-a80b-495e-886c-e51b0c7d77a3"
version = "0.1.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "d852eba0cc08181083a58d5eb9dccaec3129cb03"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.9.0"

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

    [deps.Revise.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─5a137355-8c63-492d-961c-b0fcb4897eae
# ╟─47811823-7c56-40b5-805f-8f76b0efab59
# ╟─adeb9a16-52c5-4e96-81fa-19fca745ad2b
# ╟─0ee4cf56-090d-4b9b-a266-a592952e2c1e
# ╟─14bfa5a2-aa27-4c05-a1da-d4be8037139d
# ╟─a615463a-fcaf-4823-8af4-8a04c79428ab
# ╠═0cac30f1-9101-4ba3-accc-52621bc1d16f
# ╟─0bad3ea9-d820-485e-8e84-387bf12a3f0e
# ╠═b313fdb0-3a5a-4999-a235-93cf20b3222a
# ╟─0f987ec5-7cd0-4370-b72a-2ba22e7f4c09
# ╟─3da8b6b5-f5a0-49f3-a03d-36c87584e583
# ╠═b7a2300a-bc8f-4387-a324-7f3de70d81fc
# ╠═b7880c1c-dbd2-4bf3-b916-6165c27c7d81
# ╠═d52feb78-6ddf-4b75-995e-4b2a1e452571
# ╠═8cb6660f-7181-478f-8c2e-e67c908ad8b7
# ╠═b5b5c7a1-3108-4a39-929b-f7af515b55c5
# ╠═11fa4c68-4fd8-496e-8603-8181dd3940a4
# ╠═283caef4-4cac-44fd-8fe9-17507f14baf9
# ╠═524b7699-c0f2-4904-993c-7943de80e821
# ╠═6c2ec7e6-1868-4424-ac5a-0d3ee97e0979
# ╠═760b9821-5596-46e5-b263-8101d0b2ea42
# ╠═65c9b6ea-de2f-4051-9305-6219cc814bd6
# ╠═e58e91a2-7f10-460d-9abe-155cd98e36b3
# ╠═d589347d-e0f4-40f7-b79b-f10821e19534
# ╠═b5f7928e-0299-426c-86c5-bf620bd3da03
# ╠═4347c02f-f33f-4a9c-a292-15689edafce1
# ╠═51133208-6ebf-41c4-bf52-7f5d8e294d0c
# ╠═a22422ea-683b-4d38-ac2b-d9802c0c2c80
# ╠═c1fed078-6b90-4258-a4be-87655ec3907c
# ╠═a26d7950-3b81-4a09-84ab-5e2c8fb4f00c
# ╠═0cbb83e3-7ec0-4498-9777-d6a3f56a422e
# ╠═f49155de-8db9-45b1-84f2-a85ed55a0dc5
# ╠═8bd09c6b-fedb-4256-aa10-ef6f76e9397f
# ╠═46e2e3bd-0519-4d0c-a748-fcff17118b6e
# ╟─d293006e-7bd4-46dc-82ec-dcabb7ccc5ed
# ╠═0d5207be-9594-4d28-8fa5-acd0cd397acf
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
