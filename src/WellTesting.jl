module WellTesting

using RecipesBase
using Plots
import LinearSolvers.Reg
using Statistics

export CentralDiff

include("CentralDiff.jl")
include("BuildUpRecipes.jl")
include("DrawDownRecipes.jl")




end # module
