using Revise
using LambertsProblem
using StaticArrays, LinearAlgebra

# Initialing Test Variables
r1  = SA[1.0, 0.0, 0.0]
r2  = SA[2.0*cosd(135), 2.0*sind(135), 0]
μ   = 1.0
tof = 2π

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