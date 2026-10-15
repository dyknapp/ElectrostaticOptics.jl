using ElectrostaticOptics
using StaticArrays
using Plots

# # Using `AbstractCurve2D`
# ## Line segment
# Easy, just give start and end points.

p₁ = SVector(0.0, 2.0)
p₂ = SVector(1.0, 1.0)

line = ElectrostaticOptics.LineSeg(p₁, p₂)
plot(line, aspect_ratio = :equal)

# ## CircArc
# Circular arc segment.  Give center, radius, and start/end angles

o = SVector(1.0, 2.0)
r = 1.0
θ₁ = 0.0
θ₂ = 1.75 * π

arc = CircArc(o, r, θ₁, θ₂)

plot(arc, aspect_ratio = :equal)

# Let's check the measure, should be $2\pi r$.
measure(arc, 0.0), r * abs(θ₂ - θ₁)


# ## FuncCurve
# You can define any function you like, e.g. a hyperbola.

f(t) = SVector(cosh(t), sinh(t))
hyperbola = FuncCurve(f, -1.0, 1.0)

plot(hyperbola, aspect_ratio = :equal)


# ## Test Filleting (`solve_fillet`)
# For straight lines, the fillet solve is trivial and should converge in one iteration (as long as the fillet radius is shorter than both lines).

p_corner = SVector(0.0, 0.0)

lineA = LineSeg(SVector(1.0, -1.0), p_corner)
lineB = LineSeg(p_corner, SVector(1.0, 1.0))

R = 0.25
trims, center = solve_fillet(R, (lineA, lineB))

plot(lineA, aspect_ratio = :equal)
plot!(lineB)

trim_pt_A = lineA(trims[1])
trim_pt_B = lineB(trims[2])

vec_A = trim_pt_A - center
vec_B = trim_pt_B - center

θ₁ = atan(vec_A[2], vec_A[1])
θ₂ = atan(vec_B[2], vec_B[1])

# TODO: arcs are only CCW at the moment, fix that.
fillet_arc = CircArc(center, R, θ₁, θ₂)

plot!(fillet_arc, linecolor = :blue)
