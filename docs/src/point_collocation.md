# General Overview
For a "planar kernel", we are assuming that the system extends along the $z$ axis a very long distance compared to the length scale of our system's cross-section in the $x$-$y$ plane.
We can't actually make it infinite, because then the electrostatic energy diverges, so we just take the limit.
The derivation follows [Renau1982](@cite), but our approach is a little bit different.
They start out by integrating over the full source element to find the aggregate contribution from it at a specific test (collocation) point.
This works for a collocation-specific code, but it doesn't generalize as nicely when we want to be able to take the same kernel and take a Galerkin BEM approach.
We don't want to hard-program the integral, just the point kernel.

## Derivation
We start with the classic derivation for the potential from an infinite line of charge.
Using Gauss's theorem, we know that the electric field flux leaving a cylindrical surface of height $h$ and radius $r$, centered along the line charge with charge density $\lambda$ is
```math
2\pi r h E=\frac{\lambda h}{\epsilon_0}\ .
```
The electric field strength at distance $r$ from the line charge is therefore
```math
E=\frac{\lambda}{2\pi r\epsilon_0}\ ,
```
and the field points radially.
The potential is found by integrating inward from a reference radius $r_0$, that defines the zero potential.
```math
\begin{split}
\phi(r)&=\int\vec{F}\cdot d\vec{x}\\
&=\int_{r_0}^{r}-\frac{\lambda}{2\pi\epsilon_0 R}dR\\
&=\left[\frac{\lambda}{2\pi\epsilon_0}\log\left(\frac{1}{r}\right)\right]_{r_0}^{r}\\
&=\frac{\lambda}{2\pi\epsilon_0}\log\left(\frac{r_0}{r}\right)\\
&=2\lambda k_e\log\left(\frac{r_0}{r}\right)
\end{split}
```
As you can see, taking $r_0\to\infty$ causes the potential to diverge, so we have no choice but to keep this finite reference radius hanging around.
The point kernel should give the potential contribution per unit charge, so the formula that we actually program is (K for kernel)
```math
K(r)=2k_e\log\left(\frac{r_0}{r}\right)
```
In principle, we could drop the constants out, but it's nice to know what the matrix elements of the BEM matrix physically will be.