# =====================================================================
# === Hooking into interfaces

function solve(prob::BallisticLambertsProblem{T}, ::Vallado; itermax::Int=30, tol::Float64=1e-12) where T<:Real

    @warn "Vallado solver is not fully implemented!"

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
    # longway::Bool,
    maxiter::Int,
    tol::Real
    )
	
	# Handling inputs
	rᵢ = norm(r⃗ᵢ)
	rⱼ = norm(r⃗ⱼ)

	# Finding useful constants
	δν = (retrograde ? -1 : 1)*acos( r⃗ᵢ⋅r⃗ⱼ / rᵢ / rⱼ )
    A  = (δν > π ? 1 : -1)*sqrt( (rᵢ*rⱼ)/(1+cos(δν)) )#*sin(δν)

	# Root-solving universal variable
	ψₗ, ψᵤ = -4π, 4*(π^2)
    ψ̅      = 0.0
	for k in 1:maxiter
        # Finding midpoint
        ψ̅  = 0.5*( ψₗ + ψᵤ )
        # fψ̅ = _vallado_root(ψ̅, rᵢ, rⱼ, A, μ, tof)
        c2, c3 = _c2c3(ψ̅)
        y = _y(ψ̅, c2, c3, A, rᵢ, rⱼ)

        
        if y < 0
            @warn "Lower bound too low, raising" A y
            ψₗ = (ψ̅ + 3*ψₗ)/4
            continue
        end
        x      = sqrt(y/c2)
        tof′   = (x^3 * c3 + A*sqrt(y))/sqrt(μ)

        if abs(tof′ - tof) < tol
            break
        elseif tof′ < tof
            ψₗ = ψ̅
        else
            ψᵤ = ψ̅
        end
    end
	y = _y(ψ̅, _c2c3(ψ̅)..., A, rᵢ, rⱼ)

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

function _vallado_root(ψ, rᵢ, rⱼ, A, μ, tof)
    c2, c3 = _c2c3(ψ)
    y      = _y(ψ, c2, c3, A, rᵢ, rⱼ)
    x      = sqrt(y/c2)
    tof′   = (x^3 * c3 + A*sqrt(y))/sqrt(μ)

    return tof′ - tof
end

function _y(ψ, c2, c3, A, rᵢ, rⱼ)
    # @info ψ c2 c3 A
    return rᵢ + rⱼ + ((A*(ψ*c3 - 1))/sqrt(c2))
end

function _c2c3(ψ)
    # Elliptical case
	if ψ > 1e-6
		c2 = (1.0 - cos(sqrt(ψ)))/ψ
		c3 = (sqrt(ψ) - sin(sqrt(ψ)))/sqrt(ψ^3)
		
    # Hyperbolic case
	elseif ψ < -1e-6
		c2 = (1.0 - cosh(sqrt(-ψ)))/ψ
		c3 = (sinh(sqrt(-ψ)) - sqrt(-ψ))/sqrt(-ψ^3)

    # Parabolic case
	else
		c2 = 1/2
		c3 = 1/6

	end

	return c2, c3
end

# function vallado_root(, tof)
# 	Xᵢ = √(y/c2)
# 	δtᵢ = 
# 	return δtᵢ - tof, y, δtᵢ
# end