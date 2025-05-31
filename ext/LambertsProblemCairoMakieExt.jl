module LambertsProblemCairoMakieExt

using LambertsProblem, CairoMakie, LinearAlgebra

# =====================================================================
# === Mutating plotting

function CairoMakie.lines!(ax::Axis, sol::LambertsProblem.AbstractLambertsSolution; kwargs...)
    X = getStates(sol)
    lines!(ax, [x[1] for x in X], [x[2] for x in X]; kwargs...)
end

function CairoMakie.lines!(ax::Axis3, sol::LambertsProblem.AbstractLambertsSolution; kwargs...)
    X = getStates(sol)
    lines!(ax, [x[1] for x in X], [x[2] for x in X], [x[3] for x in X]; kwargs...)
end

# =====================================================================
# === Mutating plotting

function CairoMakie.lines(sol::LambertsProblem.AbstractLambertsSolution; kwargs...)

    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())
    lines!(ax, sol; kwargs...)

    return fig, ax
end

# =====================================================================
# === Helpers

function getStates(sol::LambertsProblem.AbstractLambertsSolution)
    # Finding if hyperbolic solution$
    r, v = norm(sol.r⃗₁), norm(sol.v⃗₁)
    ε    = (0.5*v^2 - sol.μ/r)

    # Getting points
    tspan = LinRange(0.0, sol.Δt, 100 + 100*sol.revs)
    if ε < 0.0
        X = [LambertsProblem._propKepTE(sol.r⃗₁, sol.v⃗₁, t, sol.μ)[1] for t in tspan]
    else
        X = [LambertsProblem._propKepTH(sol.r⃗₁, sol.v⃗₁, t, sol.μ)[1] for t in tspan]
    end

    return X
end


end