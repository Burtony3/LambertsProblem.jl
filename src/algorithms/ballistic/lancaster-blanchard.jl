# =====================================================================
# === Hooking into interfaces

function solve(prob::BallisticLambertsProblem{T}, ::LancasterBlanchard; itermax::Int=30, tol::Float64=1e-12) where T<:Real

    # Checking inputs
    @error "This method is not implemented yet"

    # Getting solution
    # _, v⃗₁, _, v⃗₂ = _russell(prob.r⃗₁, prob.r⃗₂, prob.Δt, prob.μ, prob.retrograde, itermax, tol)

    # TODO: Calculate error
    Δtₑ = 0.0

    # Passing to output structure
    # return BallisticLambertsSolution(prob.r⃗₁, v⃗₁, prob.r⃗₂, v⃗₂, prob.Δt, Δtₑ, prob.μ, revs=prob.revs, retrograde=prob.retrograde, longway=prob.longway)

    return nothing
end


# =====================================================================
# === Main Method

#= 
# Multi-rev [Lancaster]
function _lbg(
    x⃗1::Vector{<:Real}, 
    x⃗2::Vector{<:Real}, 
    tof::Real, 
    μ::Real;
    revs::Int = 0,
    iterMax::Int = 100,
    tol::Real = 1e-12
    )
	
	# HANDLING INPUTS
	r⃗1 = x⃗ᵢ[1:3]
	v⃗1 = x⃗ᵢ[4:6]
	r1 = norm(r⃗1)
	r⃗2 = x⃗ⱼ[1:3]
	v⃗2 = x⃗ⱼ[4:6]
	r2 = norm(r⃗2)

    # Handy Constants
    r1xr2    = r⃗1×r⃗2
    r1xr2mag = norm(r1xr2)
    θ̂1       = r1xr2/r1xr2mag×r⃗1/r1
    θ̂2       = r1xr2/r1xr2mag×r⃗2/r2
    δθ       = acos(max(-1, min(1, (r⃗1⋅r⃗2)/r1/r2)))

    # Finding Branches
    islong = sign(tof)
    tof    = abs(tof)
    δθ     = islong < 0 ? δθ - 2π : δθ  
    isleft = sign(revs)
    revs   = abs(revs)

    # Defining More Constants
    c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(δθ))
    s = 0.5*(r1 + r2 + c)
    T = sqrt(8*μ/s^3)*tof
    q = sqrt(r1*r2)/s * cos(δθ/2)

end
=#