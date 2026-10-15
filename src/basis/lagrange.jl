# 2026 10 15, Daniel Knapp, Amsterdam : initial version (lagrange_j is adapted from my 2026 04 08 code)
#                                       main change is plugging into AbstractMapping


struct PolynomialMapping{T} <: AbstractMapping{1, 1}
    p::Polynomials.Polynomial{T, :x}
end
(m::PolynomialMapping)(x::Real) = m.p(x)

# 1-element vector input
function point(mapping::PolynomialMapping, ξ::SVector{1})
    return SVector(mapping.p(only(ξ)))
end
function jacobian(mapping::PolynomialMapping{T}, ξ::SVector{1}) where {T}
    deriv = derivative(mapping.p)
    return SMatrix{1, 1}(deriv(only(ξ)))
end

"""
Definition of the Lagrange basis polynomials
We construct it using the Polynomials.jl fromroots function, it's a "blink and you'll miss it" implementation.
Arguably not super readable (moving parts hidden), but it's literally the definition of the Lagrange basis.

`j` : the index of the Lagrange basis polynomial to compute (1-based indexing)
`nodes` : the vector of node positions (length = order of the polynomial)
"""
function lagrange_j(j::Int, nodes::AbstractVector{T}) where {T <: Real}
    order = length(nodes)
    @assert 1 ≤ j ≤ order "Index j must be within the bounds of the nodes"

    # We find all roots (every node except j), create the polynomial, and scale it.
    roots = [nodes[i] for i in eachindex(nodes) if i ≠ j] # Roots are all nodes except the j-th one
    p = fromroots(roots) # Construct a polynomial that has all these roots.  Still off by a constant.
    return PolynomialMapping(p / p(nodes[j])) # The polynomial is scaled to have unit height at j-th node
end
