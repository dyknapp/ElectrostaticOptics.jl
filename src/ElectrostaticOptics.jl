module ElectrostaticOptics

using ForwardDiff
using StaticArrays
using LinearAlgebra
using RecipesBase # For plot recipes

include("core/interfaces.jl")
include("geometry/abstract_curve.jl")

export AbstractCurve2D, Lineseg, CircArc, FuncCurve

end
