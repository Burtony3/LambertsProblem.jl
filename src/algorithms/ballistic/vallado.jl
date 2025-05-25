# =====================================================================
# === Hooking into interfaces

function solve(prob::BallisticLambertsProblem{T}, ::Vallado; itermax::Int=30, tol::Float64=1e-12) where T<:Real

    # Checking inputs
    @assert prob.revs == 0 "Vallado solver is a zero-revolution solver only"

    # Getting solution
    _, v⃗₁, _, v⃗₂ = _vallado(prob.r⃗₁, prob.r⃗₂, prob.Δt, prob.μ, prob.retrograde, itermax, tol)

    # TODO: Calculate error
    Δtₑ = 0.0

    # Passing to output structure
    # (r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, Δtₑ, μ, revs, retrograde, longway)
    return BallisticLambertsSolution(prob.r⃗₁, v⃗₁, prob.r⃗₂, v⃗₂, prob.Δt, Δtₑ, prob.μ, revs=prob.revs, retrograde=prob.retrograde, longway=prob.longway)
end



# =====================================================================
# === Main Method

function _vallado(
    r⃗ᵢ::SVector{3, <:Real}, 
    r⃗ⱼ::SVector{3, <:Real}, 
    tof::Real, 
    μ::Real,
    retrograde::Bool,
    maxiter::Int,
    tol::Real
    )
	
	# HANDLING INPUTS
	rᵢ = norm(r⃗ᵢ)
	rⱼ = norm(r⃗ⱼ)

	# FINDING POSITIONAL ANGLES
	δν = atand(r⃗ⱼ[2], r⃗ⱼ[1]) - atand(r⃗ᵢ[2], r⃗ᵢ[1])
	dir = retrograde ? -sign(δν) : sign(δν)

	# ETC
	cδν = r⃗ᵢ⋅r⃗ⱼ / (rᵢ*rⱼ)
	A = dir*√(rᵢ*rⱼ*(1 + cδν))
	ψₗ, ψᵤ   = -1.5, 1.5
    fψₗ, fψᵤ = vallado_root(ψₗ, rᵢ, rⱼ, A, μ, tof)[1], vallado_root(ψᵤ, rᵢ, rⱼ, A, μ, tof)[2]
    ψ̅ = 0.0
	for k in 1:maxiter
        # Finding midpoint
        ψ̅  = 0.5*( ψₗ + ψᵤ )
        fψ̅ = vallado_root(ψ̅, rᵢ, rⱼ, A, μ, tof)[1]

        # Branching
        if abs(fψ̅) < tol
            break
        elseif sign(fψ̅) == sign(fψᵤ)
            ψᵤ, fψᵤ = ψ̅, fψ̅
        elseif sign(fψ̅) == sign(fψₗ)
            ψₗ, fψₗ = ψ̅, fψ̅
        end
    end
	y, δtᵢ = vallado_root(ψ̅, rᵢ, rⱼ, A, μ, tof)[2:end]

	# f & g FUNCTIONS
	f = 1 - y/rᵢ
	g = A*√(y/μ)
	ġ = 1 - y/rⱼ

	# FINDING FINAL VALUES
	v⃗ᵢ = (r⃗ⱼ - f*r⃗ᵢ) / g
	v⃗ⱼ = (ġ*r⃗ⱼ - r⃗ᵢ) / g

	# OUTPUTTING
	return r⃗ᵢ, v⃗ᵢ, r⃗ⱼ, v⃗ⱼ
end

# =====================================================================
# === Assist Functions

function c2c3(ψ)
	if ψ > 1e-6
		c2 = (1.0 - cos(√ψ))/ψ
		c3 = (√ψ - sin(√ψ))/√(ψ^3)
		
	elseif ψ < -1e-6
		c2 = (1.0 - cosh(√(-ψ)))/ψ
		c3 = (sinh(√(-ψ)) - √(-ψ))/√(-ψ^3)

	else
		c2 = 1/2
		c3 = 1/6

	end

	return c2, c3
end

function vallado_root(ψ, rᵢ, rⱼ, A, μ, tof)
	c2, c3 = c2c3(ψ)
	y = rᵢ + rⱼ + ((A*(ψ*c3 - 1))/√(c2))
	Xᵢ = √(y/c2)
	δtᵢ = (Xᵢ^3 * c3 + A*√(y))/√(μ)
	return δtᵢ - tof, y, δtᵢ
end