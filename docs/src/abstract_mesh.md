# Overview
The "mesh" is not really the "mesh" in the sense of a 3D model.
It's not so important to know which geometry is where, but rather to be able to quickly know if two patches (that actually contain geometrical information) are connected with each other.
That allows us to dispatch quadrature for adjacent elements to a numerical integration routine that is equipped to handle the edge divergence that can ocurr between these two elements.
It's basically a list of patches, information about their connectivity, and a quick way to access the geometrical information during matrix assembly.

The basic idea for how to organize the geometry and topology information.

## Storing the actual elements
Overall: `abstract type AbstractElement{N<:Integer,D<:Integer,T<:Real} end`
Here, `N` is the dimension of the space (so 2D or 3D space), `D` is the dimension of the boundary (I think it should alost always be N-1), and `T` is just the type of real number to use.

It should support a few basic geometry functions:
`(e::AbstractElement)(\xi::T)::SVector{N,T} where {N,T}`
$\xi\in [0,1]$ is the parameter representing the coordinate in the "unit cell" that gets mapped to real geometry.
TODO is to figure out a nice way to allow multiple $\xi$, e.g. in 3D.

`jacobian(e::AbstractElement, \xi::T)::SVector{N,D,T} where {N,D,T}`
This is the Jacobian for the mapping from the unit cell to the actual geometry

The rest of the functions don't even need individual implementations for each AbstractElement: they can be derived from the Jacobbian.  Maybe nice multiple dispatch to take J directly as input (to avoid re-calculating J multiple times)
`gram(e::AbstractElement, \xi::T)::SVector{N,N,T} where {N,T}`
Gives the gram matrix `J' * J`.

`measure(e::AbstractElement, \xi::T)::SVector{N,N,T} where {N,T}`
The measure for integration.  Gives `sqrt(gram(e,\xi))`

`surface_grad(e::AbstractElement, \xi::T)::SVector{N,T} where {N,T}`
Surface gradient, for Maue regularized hypersingular kernel later on.

`normal(e::AbstractElement, \xi::T)::SVector{N,T} where {N,T}`
Normal vector.  For N=2 and D=1, `(J[2], -J[1]) / norm(J)`

## Storing the topology
`Graphs.jl` is probably nice for storing the topology.  It makes sense to have a directed graph without weights: we care about the direction of the elements for the hypersingular kernel.  `SimpleDiGraph` seems like the right option.
We should make a `Mesh` struct that keeps the vector of `AbstractEelemnt`s together with the topology graph.

## Mappings for substitutions during quadrature
Even with singularity subtraction, mappings like Telles and Duffy are nice to have to speed up quadrature convergence.
The `AbstractElement` itself is basically already a mathematical mapping...  Might make sense to have an even more general `AbstractMapping` supertype that `AbstractElement` is derived from.

## Generating test geometries
A lightweight turtle implementation would be cute/fun.