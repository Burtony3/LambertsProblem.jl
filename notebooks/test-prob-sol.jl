using Revise
using LambertsProblem
using StaticArrays, LinearAlgebra

# Initialing Test Variables
r1  = SA[1.0, 0.0, 0.0]
r2  = SA[2.0*cosd(135), 2.0*sind(135), 0]
μ   = 1.0
tof = 3π

# Creating Problem
prob = BallisticLambertsProblem(r1, r2, tof, μ, revs=0)

# Solving Problem
sol1 = solve(prob, Izzo(), tol=1e-12)
sol2 = solve(prob, Vallado(), tol=1e-12, itermax=200)
sol3 = solve(prob, Russell(), tol=1e-12)
sol3 = solve(prob, LancasterBlanchard(), tol=1e-12)

# Testing speed
using BenchmarkTools
@benchmark solve(prob, Izzo()) setup=(prob = BallisticLambertsProblem(r1, r2, 2π+randn() , μ))
@benchmark solve(prob, Vallado(), itermax=100) setup=(prob = BallisticLambertsProblem(r1, r2, 2π+randn() , μ))

# =====================================================================
# === Minimum energy problem

# Setting up problem
r1  = SA[1.0, 0.0, 0.0]
r2  = SA[2.0*cosd(165), 2.0*sind(165), 0]
μ   = 1.0
prob = MinimumEnergyLambertsProblem(r1, r2, μ, retrograde=false)

# Solving
sol = solve(prob)

# Testing against izzo
prob = BallisticLambertsProblem(r1, r2, sol.Δt, μ, revs=sol.revs, retrograde=sol.retrograde, longway=false)
sol2 = solve(prob, Izzo(), tol=1e-14)

sol.v⃗₁ - sol2.v⃗₁
sol.v⃗₂ - sol2.v⃗₂

# =====================================================================
# === Plotting

using CairoMakie

# Getting values to plot
r1  = SA[1.0, 0.0, 0.0]
θ   = 135.0
r2  = SA[2.0*cosd(θ), 2.0*sind(θ), 0]
μ   = 1.0
tof = 0.5π

# Solving minimum energy problem
prob = MinimumEnergyLambertsProblem(r1, r2, μ)
sol_le  = solve(prob)
prob = MinimumEnergyLambertsProblem(r1, r2, μ, revs=1)
sol_ler  = solve(prob)

# Solving fixed-tof lambert
prob = BallisticLambertsProblem(r1, r2, tof, μ, longway=false)
sol_b = solve(prob, Izzo())
prob = BallisticLambertsProblem(r1, r2, 8π, μ, revs=1, longway=false)
sol_br = solve(prob, Izzo())
prob = BallisticLambertsProblem(r1, r2, 8π, μ, revs=1, longway=true)
sol_brl = solve(prob, Izzo())
prob = BallisticLambertsProblem(r1, r2, sol_le.Δt, μ)
sol_ble = solve(prob, Izzo())

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

# Plotting
begin
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    # Plotting 
    plot_sol!(ax, sol_le)
    # plot_sol!(ax, sol_ler)
    plot_sol!(ax, sol_b)
    plot_sol!(ax, sol_br)
    plot_sol!(ax, sol_brl)
    plot_sol!(ax, sol_ble, linestyle=:dash)

    fig
end