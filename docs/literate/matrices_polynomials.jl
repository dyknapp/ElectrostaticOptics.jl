# TODO: this is a stream-of-thought documentation file.  Rewrite nicely.

using LinearAlgebra
using Polynomials
using ElectrostaticOptics
using Plots

# # General idea
# A polynomial is just a vector in the polynomial vector space.  For a vector $\vec{c}$:
# ```math
#   P(\vec{c}, x)=c_1+c_2x+c_3x^2+c_4x^3+\dots
# ```
# We can do basis changes to use different polynomial bases.  E.g., in the Lagrange polynomial basis, we can take a vector $\vec{l}$
# and go back to the naive "coefficient space" vector $\vec{c}$ as
# ```math
#   \vec{c}=\mathbf{L}\vec{l}\ .
# ```
# A derivative is just a linear operator, represented by the "coefficient space" matrix,
# TODO: fix this formula
# ```math
#   \mathbf{D}=
#   \begin{pmatrix}
#       \ddots & 1 &   &   &   \\
#       \ddots & 0 & 2 &   &   \\
#              &   & 0 & 3 &   \\
#              &   &   & 0 &   \\
#              &   &   &   & 0
#   \end{pmatrix}\ .
# ```
function derivative_matrix(n::Int)
    D = zeros(Float64, n, n)
    for idx in 2:n
        D[idx - 1, idx] = Float64(idx - 1)
    end

    return D
end

D = derivative_matrix(4)

## Let's test it out

p = Polynomial([1, 2, 3, 4])

## Derivative should be $12x^2+6x+2$

Polynomial(derivative_matrix(4) * p.coeffs)

## We can also look at it in the Lagrange basis

function lagrange_basis_matrix(n::Int)
    L = zeros(Float64, n, n)
    lb = LagrangeBasis(Float64, n - 1)
    for (col, poly) in enumerate(lb.mappings)
        L[:, col] = poly.p.coeffs
    end

    return L
end

L = lagrange_basis_matrix(4)

## Great, let's try taking a derivative, I'm going to make up a random polynomial in the Lagrange basis

l = [1.0, 0.0, 3.0, 0.0]

lp = Polynomial(L * l)

## The derivative should be $-19+126x-135x^2$.

Polynomial(D * L * l)

## We can just do a basis change if we want to find the derivatives in another basis:
# ```math
#   \vec{dl}=\mathbf{L}^{-1}\mathbf{D}\mathbf{L}\vec{l}\ .
# ```

dlp = (inv(L) * D * L) * l

## We can easily check that this matches the expected derivative (within floating-point error):

Polynomial(L * dlp)


## Evaluation of a polynomial takes place using a vector of powers $\vec{p}=[x^0, x^1, x^2, x^3, \dots]$

function power_vector(n::Int, x::Float64)
    return transpose((x) .^ collect(0:(n - 1)))
end

power_vector(4, 3.0)

## Let's try it out with a simple polynomial.  We expect:

ptest = Polynomial([1.0, 2.0, 3.0, 4.0])
ptest(3.14)

## We get the same with the power vector:

power_vector(4, 3.14) * ptest.coeffs

## In general, if we want a vector of the polynomial evaluated at $n$ points, we'd get the Vandermonde matrix, leading to the polynomial interpoaltion problem.
# Anyways, what we really want is to constrain our DOFs so by imposing slope or value conditions at the endpoints.
# Let's say for instance that we are using cubic Lagrange elements, and we want to clamp the endpoints at specific values.
# This is pretty trivial, but we want to see it falling out of the linear algebnra.
# The endpoint clamps give us two conditions, but we need 4 to determine the coefficients, so we need to keep two internal DOFs.
# For simplicity, let's put those two DOFs at $\xi=1/4,3/4$, just for simplicity (it's easy to see on a plot).
# Ok, let's derive our system of equations.

D = derivative_matrix(4)
L = lagrange_basis_matrix(4)
ξ0 = power_vector(4, 0.0) # all zeros except for the first element (0^0=1)
ξ1_4 = power_vector(4, 1.0 / 4.0)
ξ3_4 = power_vector(4, 3.0 / 4.0)
ξ1 = power_vector(4, 1.0) # all ones

# Row: derivative at ξ=0
r1 = ξ0 * (D * L)

# Row: derivative at ξ=1
r2 = ξ1 * (D * L)

# Row: value at ξ=1/3
r3 = ξ1_4 * L

# Row: value at ξ=2/3
r4 = ξ3_4 * L

# The full matrix
M = [r1; r2; r3; r4]

## Okay, let's find a polynomial with $10.0$ slope at the end points and the mid points at $\mp4.0$

# Solve for Lagrange coefficients
l_sol = M \ [10.0, 10.0, 4.0, -4.0]

# Map Lagrange coefficients to standard polynomial coefficients: c = L * l
result_polynomial = Polynomial(L * l_sol)

# Verify constraints
dp = derivative(result_polynomial)
println("p'(0)   = ", dp(0.0))
println("p'(1)   = ", dp(1.0))
println("p(1/4)  = ", result_polynomial(1 / 4))
println("p(3/4)  = ", result_polynomial(3 / 4))

x = LinRange(0.0, 1.0, 128)
plot(x, result_polynomial.(x), label = "p(x)", legend = false)

# In general, this is how we can turn a chosen polynomial basis into global DOFs.
# At each joint between panels, we can choose a set number of constraints (e.g., for smooth charge distribution, two constraints, derivative and value).
# The DOF assignment then proceeds by allocating unconstrained DOFs to the interio  r of each panel, or erroring if polynomial order is too low to satisfy constraints.
# The topology graph edges request constraints, and the vertices soak up the remaining ones.
# This gives an easy way to give raw basis functions directly corresponding to global DOFs to the panels.
