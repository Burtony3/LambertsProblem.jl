# =====================================================================
# === Hooking into interfaces

function solve(prob::BallisticLambertsProblem{T}, ::Russell; itermax::Int=30, tol::Float64=1e-12) where T<:Real

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