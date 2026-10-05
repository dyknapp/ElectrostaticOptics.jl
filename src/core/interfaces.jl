# 2026 10 04, Daniel Knapp, Amsterdam : initial version (mostly writing out conventions)


# AbstractSpace
#
# A "space" (TODO: think of / find a better name) is the choice of 2 coordinates to use as our 2D problem setup.
# For instance, planar problems will use (x, y) and axisymmetric problems will use (r, z).
# The idea behind making this an abstract supertype is that I want to be able to use various types of patches (but not kernels!) in different spaces.
# The problem I had before is that the Jacobian was hard-coded into the patch, which was irritating when I wanted to switch from axisymmetric to planar.
# What I'm currently wondering is what I actually want to build into AbstractSpace.  I can imagine some nice and simple stuff like definitions for distances or whatever.
# But is it going to just end up being a label for "which Jacobian should I use?".  Maybe.
# Not an important design decision immediately because I'm initially going to focus on point collocation BEM.  But something to keep in mind!!
#
abstract type AbstractSpace end

# AbstractPatch
#
# A patch is a 1D curve in 2D space, defined by a geometry map from ξ∈[0,1] to 2D coordinates (e.g. (x,y) for planar or (r,z) for axisymmetric).
# As a basic design convention, we follow a hierarchical structure: each surface charge basis function must have support on a single patch.
#   EXCEPTION: adjacent patches can share a basis function that has support over both, as a way to enforce continuity.
#       The basic rule is one DOF per control point.  Control points belong to patches but can be shared.
# A lot of the basis function computation logic is attached to the patch because of this heirarchy:
#   A lot of information (polynomial coefficients, geometrical information like the Jacobian) is inherently tied to the geometric patch
#
#
# Abstract interface for patches:
#
# - patch(ξ) → SVector{2}  (evaluate coordinate corresponding to ξ)
#       This one is self-explanatory, we need to know the geometry...
#
# - jacobian(patch, ξ) → Real (evaluate the Jacobian at ξ)
#       Also relatively self-explanatory, we need the Jacobian for integration.
#
# - control_points(patch) → Vector{SVector{2}}  (return the control points or nodes of the patch)
#       This is needed for plotting and for adjacency detection.
#           By default, the geometry is loaded into this package without any information on adjacency or ordering:
#           During integration, we need to know adjacency to to predict when an integrand will be singular and need special treatment.
#
abstract type AbstractPatch{T <: Real, N} end
(patch::AbstractPatch)(ξ::Real) = error("Not implemented for $(typeof(patch))")
jacobian(patch::AbstractPatch, space::AbstractSpace, ξ::Real) = error("Not implemented for $(typeof(patch))")
control_points(patch::AbstractPatch) = error("Not implemented for $(typeof(patch))")

# AbstractKernel
#
# A kernel defines the interaction between two points in space, e.g. the Coulomb interaction for electrostatics.
# Since this package has a focus on 2D BEM, even for electrostatics (or other possible future things), the kernels will look a bit funny
#   because we're going to be integrating out various degrees of freedom.
#
# Abstract interface for kernels:
#
# - patch(pt₁, pt₂) → Real (evaluate the kernel between two points)
#       Just evaluate the formula between two points.
#       WARNING: this will almost certainly be infinite if pt₁ = pt₂.
#           You probably want to avoid that happening in your code (singular quadrature), or handle it explicitly (limiting value for point collocation)
abstract type AbstractKernel{S <: AbstractSpace} end
(kernel::AbstractKernel, pt₁, pt₂) = error("Not implemented for $(typeof(kernel))")
# TODO: is there a nice general way to ask for diagonal/adjacent element/far separated element blocks?
