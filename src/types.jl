# =====================================================================
# === Problems

abstract type AbstractLambertsProblem end

# Minimum Energy Problem
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

# Ballistic Problems
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

abstract type AbstractLambertsSolvers end

struct Izzo <: AbstractLambertsSolvers end
struct Vallado <: AbstractLambertsSolvers end
struct Russell <: AbstractLambertsSolvers end
struct LancasterBlanchard <: AbstractLambertsSolvers end

export AbstractLambertsSolvers,
       Izzo,
       Vallado,
       Russell,
       LancasterBlanchard


# =====================================================================
# === Solutions

abstract type AbstractLambertsSolution end

# Ballistic Solution
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

# Ballistic Solution
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