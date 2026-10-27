# 2026 10 15, Daniel Knapp, Amsterdam : initial version (lagrange_j is adapted from my 2026 04 08 code)
#                                       main change is plugging into AbstractMapping
# 2026 10 26, Daniel Knapp, Amsterdam : LagrangeBasis constructor for PolynomialBasis

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

# https://discourse.julialang.org/t/test-whether-an-iterable-is-strictly-increasing/100282
function is_pairwise(f, itr)
    p = Iterators.peel(itr)
    isnothing(p) && return true
    prev, rest = p
    for cur in rest
        f(prev, cur) || return false
        prev = cur
    end
    return true
end

# e.g. for P=2, quadratic so 3 degrees of freedom
function LagrangeBasis(local_nodes::SVector{N, T}) where {N, T <: Real}
    @assert N ≥ 1
    @assert all(x -> 0 ≤ x ≤ 1, local_nodes)
    @assert is_pairwise(<, local_nodes)

    mappings = SVector{N, PolynomialMapping}(
        ntuple(j -> lagrange_j(j, local_nodes), N)
    )

    return PolynomialBasis{N, PolynomialMapping{T}}(mappings)
end

function LagrangeBasis(::Type{T}, p::Int) where {T <: Real}
    @assert p ≥ 0
    if p == 0
        local_nodes = [0.5]
    else
        local_nodes = LinRange(T(0), T(1), p + 1)
    end
    mappings = SVector{p + 1, PolynomialMapping}(
        ntuple(j -> lagrange_j(j, local_nodes), p + 1)
    )
    return PolynomialBasis{p + 1, PolynomialMapping{T}}(mappings)
end

# Default: Float64
LagrangeBasis(p::Int) = LagrangeBasis(Float64, p)


@recipe function f(pb::PolynomialBasis)
    seriestype   --> :path
    labels       --> false
    ylims        --> (-1.5, 1.5)

    for idx in eachindex(pb.mappings)
        @series begin
            t = LinRange(0.0, 1.0, 128)
            pts = pb.mappings[idx].(t)
            t, pts
        end
    end
end
