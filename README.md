# LambertsProblem

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://burtony3.github.io/LambertsProblem.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://burtony3.github.io/LambertsProblem.jl/dev/)
[![Build Status](https://github.com/burtony3/LambertsProblem.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/burtony3/LambertsProblem.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/burtony3/LambertsProblem.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/burtony3/LambertsProblem.jl)

A Julia package providing a variety of algorithms to set up and solve the orbital mechanics Lambert problem.

---

## Features

* Define ballistic Lambert problems in ℝ³ with support for multiple revolutions, retrograde, and/or short/long transfers.
* Multiple solver implementations:

  * Izzo
  * Vallado (incomplete)
  * Russell (not implemented)
  * Lancaster–Blanchard (not implemented)

* Tolerance and iteration control
* Support for multiple revolutions

## Installation

```julia
using Pkg
Pkg.add("LambertsProblem")
```

## Usage

```julia
using LambertsProblem
using StaticArrays, LinearAlgebra

# Initialize test vectors
r⃗₁ = SA[1.0, 0.0, 0.0]
r⃗₂ = SA[2.0*cosd(135), 2.0*sind(135), 0.0]
μ   = 1.0                 # Standard gravitational parameter
Δt  = 3π                  # Time-of-flight

# Create a ballistic Lambert problem
prob = BallisticLambertsProblem(
    r⃗₁, r⃗₂, Δt, μ;
    revs=0, retrograde=false, longway=false
)

# Solve with different algorithms
sol_izzo      = solve(prob, Izzo(); tol=1e-12)
sol_vallado   = solve(prob, Vallado(); tol=1e-12, itermax=200)
sol_russell   = solve(prob, Russell(); tol=1e-12)
sol_lancaster = solve(prob, LancasterBlanchard(); tol=1e-12)
```

Each solution `sol` contains the velocity vectors at departure and arrival, plus metadata:

```julia
sol.v⃗₁           # departure velocity
sol.v⃗₂           # arrival velocity
sol.Δt           # solved time-of-flight
sol.Δtₑ          # time-of-flight error from requested
```

## Solvers

* **Izzo**: Fast general-purpose, multi-rev solver
* **Vallado**: Classical zero-rev solver
* **Russell**: Alternative multi-rev solver
* **Lancaster–Blanchard**: Robust solver

## API Reference

### Core Types

```julia
abstract type AbstractLambertsProblem end

struct BallisticLambertsProblem{T<:Real} <: AbstractLambertsProblem
    r⃗₁::SVector{3,T}
    r⃗₂::SVector{3,T}
    Δt::T
    μ::T
    revs::Int
    retrograde::Bool
    longway::Bool
    
    function BallisticLambertsProblem(
        r⃗₁::SVector{3,T}, r⃗₂::SVector{3,T},
        Δt::T, μ::T;
        revs::Int=0, retrograde::Bool=false, longway::Bool=false
    ) where T<:Real
        @assert Δt>0 "Time of flight must be positive"
        @assert μ>0  "Standard gravitational parameter must be positive"
        @assert revs≥0 "Number of revs must not be negative"
        return new{T}(r⃗₁, r⃗₂, Δt, μ, revs, retrograde, longway)
    end
end
```

```julia
abstract type AbstractLambertsSolvers end
struct Izzo             <: AbstractLambertsSolvers end
struct Vallado          <: AbstractLambertsSolvers end
struct Russell          <: AbstractLambertsSolvers end
struct LancasterBlanchard <: AbstractLambertsSolvers end
```

```julia
abstract type AbstractLambertsSolution end

struct BallisticLambertsSolution{T<:Real} <: AbstractLambertsSolution
    r⃗₁::SVector{3,T}
    v⃗₁::SVector{3,T}
    r⃗₂::SVector{3,T}
    v⃗₂::SVector{3,T}
    Δt::T
    Δtₑ::T
    μ::T
    revs::Int
    retrograde::Bool
    longway::Bool

    function BallisticLambertsSolution(
        r⃗₁::SVector{3,T}, v⃗₁::SVector{3,T},
        r⃗₂::SVector{3,T}, v⃗₂::SVector{3,T},
        Δt::T, Δtₑ::T, μ::T;
        revs::Int=0, retrograde::Bool=false, longway::Bool=false
    ) where T<:Real
        @assert Δt>0 "Time of flight must be positive"
        @assert μ>0  "Standard gravitational parameter must be positive"
        @assert revs≥0 "Number of revs must not be negative"
        return new{T}(r⃗₁, v⃗₁, r⃗₂, v⃗₂, Δt, Δtₑ, μ, revs, retrograde, longway)
    end
end
```

### Solver Interface

```julia
solve(
    prob::BallisticLambertsProblem,
    algo::AbstractLambertsSolvers;
    tol::Real=1e-10,
    itermax::Int=1000
) -> BallisticLambertsSolution
```

* `tol`: convergence tolerance.
* `itermax`: maximum iterations (only for algorithms that support it).
* Problem-level options (`revs`, `retrograde`, `longway`) set on construction.

---

## Contributing

1. Fork the repository
2. Create a feature branch
3. Write tests under `test/`
4. Document additions
5. Open a pull request
