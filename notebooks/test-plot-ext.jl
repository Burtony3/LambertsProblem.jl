using Revise
using LambertsProblem
using StaticArrays, LinearAlgebra

# Initialing Test Variables
r1  = SA[1.0, 0.0, 0.0]
r2  = SA[2.0*cosd(135), 2.0*sind(135), 0.5]
μ   = 1.0
tof = 3π

# Creating Problem
prob = BallisticLambertsProblem(r1, r2, tof, μ, revs=0)
sol1 = solve(prob, Izzo())
prob = MinimumEnergyLambertsProblem(r1, r2, μ, revs=0)
sol2 = solve(prob)

# Plotting
using CairoMakie
fig, ax = lines(sol1)
begin
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    lines!(ax, sol1)
    lines!(ax, sol2)
 
    fig
end

begin
    fig = Figure()
    ax  = Axis3(fig[1, 1], aspect=:equal)

    lines!(ax, sol1) 
    lines!(ax, sol2) 

    fig
end