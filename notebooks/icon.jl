using Revise
using LambertsProblem
using StaticArrays, LinearAlgebra
using CairoMakie, Colors

# Initialing Test Variables
R, θ = 1.0, -25.0
r1   = SA[R*cosd(θ), R*sind(θ), 0]
R, θ = 2.0, 120.0
r2   = SA[R*cosd(θ), R*sind(θ), 0]
μ    = 1.0

# Finding minimum energy transfer
prob = MinimumEnergyLambertsProblem(r1, r2, μ)
sol_le  = solve(prob)

# Creating plotting function
function plot_sol!(ax, sol; kwargs...)
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

    # Plotting
    lines!(ax, [x[1] for x in X], [x[2] for x in X]; kwargs...)
end

# Making plot
logocolors = Colors.JULIA_LOGO_COLORS
begin
    fig = Figure(backgroundcolor = (:white, 0.0))
    ax  = Axis(fig[1, 1], aspect=DataAspect(), backgroundcolor=(:white, 0.0))
    hidedecorations!(ax)
    hidespines!(ax)
    
    # Plotting low energy transfer
    plot_sol!(ax, sol_le; color=logocolors.blue, linewidth=3)

    # Plotting other transfers
    for (dt, c) in zip([0.5π, 3π, 4.5π], [logocolors.red, logocolors.green, logocolors.purple])
        prob = BallisticLambertsProblem(r1, r2, dt, μ)
        sol  = solve(prob, Izzo())

        plot_sol!(ax, sol; color=c, linewidth=3)
    end

    # Adding other accents1
    scatter!(ax, [r1[1], r2[1]], [r1[2], r2[2]], color=:black, markersize=12)
    R = 0.25
    θ⃗ = LinRange(0.0, 2π, 100)
    poly!(ax, 
        Point2f[(R*cos(θ), R*sin(θ)) for θ in θ⃗],
        color = :dimgray, strokecolor = :black, strokewidth=3
    )
    # Label(fig[1, 2:10], "LambertsProblem.jl", fontsize=10, color=:white, 
        # justification=:left, width=3)

    fig
end

save(@__DIR__()*"/../icon.svg", fig)