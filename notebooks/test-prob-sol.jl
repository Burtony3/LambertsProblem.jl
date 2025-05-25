using Revise
using LambertsProblem
using StaticArrays, LinearAlgebra

# Initialing Test Variables
r1  = SA[1.0, 0.0, 0.0]
r2  = SA[2.0*cosd(135), 2.0*sind(135), 0]
μ   = 1.0
tof = 10π

# Creating Problem
prob = BallisticLambertsProblem(r1, r2, tof, μ, revs=1)

# Solving Problem
sol = solve(prob, Izzo())

# Testing speed
using BenchmarkTools
@benchmark solve(prob, Izzo()) setup=(prob = BallisticLambertsProblem(r1, r2, 2π+randn() , μ))