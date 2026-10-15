module ElectrostaticOptics

using ForwardDiff
using StaticArrays
using LinearAlgebra
using RecipesBase # For plot recipes

include("core/interfaces.jl")
include("geometry/abstract_curve.jl")
include("geometry/fillet.jl")

# geometry
export AbstractCurve2D, LineSeg, CircArc, FuncCurve
export tangent, normal, measure
export solve_fillet, solve_fillet_arc

end
