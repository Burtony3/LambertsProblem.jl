# =====================================================================
# === Problems

"""
    AbstractLambertsProblem

Abstract supertype for all Lambert’s‐problem definitions.
"""
abstract type AbstractLambertsProblem end

"""
    MinimumEnergyLambertsProblem(r⃗₁, r⃗₂, μ; revs=0, retrograde=false, longway=false)

Defines a minimum‐energy Lambert problem between two 3D position vectors under gravitational parameter `μ`.

- `revs`: number of allowed revolutions (default 0).  
- `retrograde=true`: force a retrograde (opposite‐direction) trajectory.  
- `longway=true`: force the longer‐path geometry, even if it costs more energy.
"""
struct MinimumEnergyLambertsProblem{T<:Real} <: AbstractLambertsProblem
    # Required arguments
    r⃗₁::SVector{3, T}
    r⃗₂::SVector{3, T}
    μ::T

    # Optional Arguements
    revs::Int
    retrograde::Bool
    longway::Bool

    # Constructor
    function MinimumEnergyLambertsProblem(
        r⃗₁::SVector{3, T}, r⃗₂::SVector{3, T}, μ::T; 
        revs::Int=0, retrograde::Bool=false, longway::Bool=false
        ) where T<:Real
        @assert μ>0 "Standard gravitational parameter must be positive"
        @assert revs≥0 "Number of revs must not be negative"
        return new{T}(r⃗₁, r⃗₂, μ, revs, retrograde, longway)
    end
end

"""
    BallisticLambertsProblem(r⃗₁, r⃗₂, Δt, μ; revs=0, retrograde=false, longway=false)

Defines a fixed‐time (ballistic) Lambert problem between two 3D vectors, solving for departure and arrival velocities that satisfy the given time‐of‐flight `Δt`.

- `revs`: number of allowed revolutions (default 0).  
- `retrograde=true`: force a retrograde (opposite‐direction) trajectory.  
- `longway=true`: force the longer‐path geometry, even if it requires a longer angle.
"""
struct BallisticLambertsProblem{T<:Real} <: AbstractLambertsProblem
    # Required arguments
    r⃗₁::SVector{3, T}
    r⃗₂::SVector{3, T}
    Δt::T
    μ::T

    # Optional Arguements
    revs::Int
    retrograde::Bool
    longway::Bool

    # Constructor
    function BallisticLambertsProblem(
        r⃗₁::SVector{3, T}, r⃗₂::SVector{3, T},
        Δt::T, μ::T; 
        revs::Int=0, retrograde::Bool=false, longway::Bool=false
        ) where T<:Real
        @assert Δt>0 "Time of flight must be positive"
        @assert μ>0 "Standard gravitational parameter must be positive"
        @assert revs≥0 "Number of revs must not be negative"
        return new{T}(r⃗₁, r⃗₂, Δt, μ, revs, retrograde, longway)
    end
end

# Impulsive Problems
# struct ImpulsiveLambertsProblem{T<:Real} <: AbstractLambertsProblem
# end

export AbstractLambertsProblem,
       MinimumEnergyLambertsProblem,
       BallisticLambertsProblem


# =====================================================================
# === Solvers

"""
    AbstractLambertsSolvers

Abstract supertype for all Lambert‐problem solver algorithms.
"""
abstract type AbstractLambertsSolvers end

"""
    Izzo()

Izzo’s universal‐variable Lambert solver (Householder iterations).
"""
struct Izzo <: AbstractLambertsSolvers end

"""
    Vallado()

Vallado’s series‐expansion Lambert solver.
"""
struct Vallado <: AbstractLambertsSolvers end

"""
    Russell()

Russell’s successive‐substitution Lambert solver.
"""
struct Russell <: AbstractLambertsSolvers end

"""
    LancasterBlanchard()

Lancaster–Blanchard analytical Lambert solver (multiple‐revolution cases).
"""
struct LancasterBlanchard <: AbstractLambertsSolvers end

export AbstractLambertsSolvers,
       Izzo,
       Vallado,
       Russell,
       LancasterBlanchard


# =====================================================================
# === Solutions

"""
    AbstractLambertsSolution

Abstract supertype for all Lambert‐problem solution containers.
"""
abstract type AbstractLambertsSolution end

"""
    MinimumEnergyLambertsSolution(r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, μ; revs=0, retrograde=false, longway=false)

Solution container for a minimum‐energy Lambert problem. Includes initial/final positions `r⃗₁`, `r⃗₂`, computed velocities `v⃗₁`, `v⃗₂`, time of flight `Δt`, and gravitational parameter `μ`.

- `revs`: number of revolutions used.  
- `retrograde=true`: indicates the solution is retrograde (opposite‐direction).  
- `longway=true`: indicates the solution follows the longer‐path geometry.
"""
struct MinimumEnergyLambertsSolution{T<:Real} <: AbstractLambertsSolution 
    # Required arguments
    r⃗₁::SVector{3, T}
    v⃗₁::SVector{3, T}
    r⃗₂::SVector{3, T}
    v⃗₂::SVector{3, T}
    Δt::T
    μ::T

    # Optional Arguements
    revs::Int
    retrograde::Bool
    longway::Bool

    # Constructor
    function MinimumEnergyLambertsSolution(
        r⃗₁::SVector{3, T}, v⃗₁::SVector{3, T},
        r⃗₂::SVector{3, T}, v⃗₂::SVector{3, T},
        Δt::T, μ::T; 
        revs::Int=0, retrograde::Bool=false, longway::Bool=false
        ) where T<:Real
        @assert Δt>0 "Time of flight must be positive"
        @assert μ>0 "Standard gravitational parameter must be positive"
        @assert revs≥0 "Number of revs must not be negative"
        return new{T}(r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, μ, revs, retrograde, longway)
    end
end

"""
    BallisticLambertsSolution(r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, Δtₑ, μ; revs=0, retrograde=false, longway=false)

Solution container for a ballistic (fixed‐time) Lambert problem. Includes computed departure/arrival velocities `v⃗₁`, `v⃗₂`, nominal time of flight `Δt`, optional ephemeris deviation `Δtₑ`, and gravitational parameter `μ`.

- `revs`: number of revolutions used.  
- `retrograde=true`: indicates the solution is retrograde (opposite‐direction).  
- `longway=true`: indicates the solution follows the longer‐path geometry.
"""
struct BallisticLambertsSolution{T<:Real} <: AbstractLambertsSolution 
    # Required arguments
    r⃗₁::SVector{3, T}
    v⃗₁::SVector{3, T}
    r⃗₂::SVector{3, T}
    v⃗₂::SVector{3, T}
    Δt::T
    Δtₑ::T
    μ::T

    # Optional Arguements
    revs::Int
    retrograde::Bool
    longway::Bool

    # Constructor
    function BallisticLambertsSolution(
        r⃗₁::SVector{3, T}, v⃗₁::SVector{3, T},
        r⃗₂::SVector{3, T}, v⃗₂::SVector{3, T},
        Δt::T, Δtₑ::T, μ::T; 
        revs::Int=0, retrograde::Bool=false, longway::Bool=false
        ) where T<:Real
        @assert Δt>0 "Time of flight must be positive"
        @assert μ>0 "Standard gravitational parameter must be positive"
        @assert revs≥0 "Number of revs must not be negative"
        return new{T}(r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, Δtₑ, μ, revs, retrograde, longway)
    end
end

# Impulsive Solution
# struct ImpulsiveLambertsSolution{T<:Real} <: AbstractLambertsSolution
# end

export AbstractLambertsSolution,
       MinimumEnergyLambertsSolution,
       BallisticLambertsSolution