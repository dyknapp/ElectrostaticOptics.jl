module ElectrostaticOptics

using ForwardDiff
using StaticArrays
using LinearAlgebra
using RecipesBase # For plot recipes
using Graphs # For mesh topology
using Polynomials

include("core/interfaces.jl")
include("geometry/abstract_curve.jl")
include("geometry/fillet.jl")
include("geometry/mesh.jl")
include("basis/lagrange.jl")
include("basis/dofs.jl")

# geometry
export AbstractCurve2D, LineSeg, CircArc, EllipticArc, FuncCurve
export tangent, normal, measure, derivative_curve
export affine_transform, translation
export solve_fillet, solve_fillet_arc

# mesh
export empty_mesh, concat_meshes, circle, rectangle
export fillet_corners

# basis
export lagrange_j, jacobian
export LagrangeBasis

export BasisConstraint, UNCONSTRAINED, VALUE, C0, C1
export build_function_space, local_coefficients

end
