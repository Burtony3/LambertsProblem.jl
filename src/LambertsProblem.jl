module LambertsProblem

# Package Imports
using LinearAlgebra, StaticArrays

# Definitions
include("types.jl")

# Solvers
include("algorithms/ballistic/izzo.jl")
include("algorithms/ballistic/vallado.jl")
include("algorithms/ballistic/russell.jl")
include("algorithms/ballistic/lancaster-blanchard.jl")

end
