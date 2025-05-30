module LambertsProblem

# Package Imports
using LinearAlgebra, StaticArrays

# Definitions
include("types.jl")
include("utils.jl")

# Solvers
include("algorithms/ballistic/min-energy.jl")
include("algorithms/ballistic/izzo.jl")
include("algorithms/ballistic/vallado.jl")
include("algorithms/ballistic/russell.jl")
include("algorithms/ballistic/lancaster-blanchard.jl")

end
