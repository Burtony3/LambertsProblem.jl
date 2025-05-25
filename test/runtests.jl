using LambertsProblem
using Test
using StaticArrays, LinearAlgebra

@testset "Ballistic Lambert Solvers — Zero Rev" begin
    # Input variables
    r1  = SA[1.0, 0.0, 0.0]
    r2  = SA[2.0*cosd(135), 2.0*sind(135), 0]
    μ   = 1.0
    tof = 2π
    prob = BallisticLambertsProblem(r1, r2, tof, μ)
    tol     = 1e-12

    # Test values
    v1′, v2′ = SA[0.3793321336652016, 1.0832820655755302, 0.0], SA[-0.2734127087511391, -0.49258338575508864, 0.0]

    # Testing lambert solvers against one another
    sol = solve(prob, Izzo(), itermax=20, tol=tol)
    @test sol.v⃗₁ ≈ v1′ atol=1e-10
    @test sol.v⃗₂ ≈ v2′ atol=1e-10
    # sol = solve(prob, Vallado(), itermax=300, tol=tol)
    # @test sol.v⃗₁ ≈ v1′ atol=1e-10 skip=true
    # @test sol.v⃗₂ ≈ v2′ atol=1e-10 skip=true
end

@testset "Ballistic Lambert Solvers — Multi Rev" begin
    # Input variables
    r1  = SA[1.0, 0.0, 0.0]
    r2  = SA[2.0*cosd(135), 2.0*sind(135), 0]
    μ   = 1.0
    tof = 10π
    prob = BallisticLambertsProblem(r1, r2, tof, μ, revs=1)
    itermax = 20
    tol     = 1e-12

    # Test values
    v1′, v2′ = SA[-0.22919931718118566, 1.2582313092559003, 0.0], SA[-0.7911840460500928, -0.09851984502598232, 0.0]

    # Testing lambert solvers against one another
    sol = solve(prob, Izzo(), itermax=20, tol=tol)
    @test sol.v⃗₁ ≈ v1′ atol=1e-10
    @test sol.v⃗₂ ≈ v2′ atol=1e-10
    @test_throws AssertionError solve(prob, Vallado(), itermax=300, tol=tol)
end
