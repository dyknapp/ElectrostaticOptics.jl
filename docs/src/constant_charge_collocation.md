# Overview
The simplest version of the BEM electrostatics solver in this package is a point collocation solver with straight line constant charge density basis elements.
This is a comparison point to benchmark the "higher-order" methods that we will be implementing and mainly relying on.

## Planar kernel
Renau, Read, and Brunt [Renau1982](@cite) analytically integrate (and approximate the result where necessary) to obtain a formula that directly gives their matrix elements.
By contrast, we will take the path of singular quadrature.
It might be a bit slower, but it will converge to the exact value (better than the approximate one that RRH give), and singular quadrature is the main challenge of higher-order methods anyways.
The planar kernel is an especially nice test case, because of its simple form, with an explicit weak (logarithmic) divergence.
In 2D BEM electrostatics, logarithmic divergences are going to be the main/only form of divergence that we see (at least until we get to hypersingular kernels for Calderon preconditioning), but it already get's more complicated with axisymmetric kernels and the elliptic integrals that come into play.

### Quadrature method
We will see logarithmic divergences for the diagonal matrix elements (the "self-interaction") and for matrix elements that correspond to adjacent elements, because one element's edge will "touch" the the edge of its neighbor.
We need in general to handle those two cases separately, so there needs to be some awareness of what the mesh is like to dispatch those matrix elements to the correct integration function.
We will have multiple logarithmic-weighted quadrature routines to benchmark against each other.

- log-weighted Gaussian quadrature
- log-weighted Clenshaw-Curtis quadrature
- tanh-sinh quadrature

