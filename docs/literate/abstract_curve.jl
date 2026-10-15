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
θ₁ = 0.5 * π
θ₂ = 0.2 * π

arc = CircArc(o, r, θ₁, θ₂)

plot(arc, aspect_ratio = :equal)

# Let's check the measure, should be $2\pi r$.
# measure(arc, 0.0), r * abs(θ₂ - θ₁)


# ## FuncCurve
# You can define any function you like, e.g. a hyperbola.

f(t) = SVector(cosh(t), sinh(t))
hyperbola = FuncCurve(f, -1.0, 1.0)

plot(hyperbola, aspect_ratio = :equal)


# ## Test Filleting (`solve_fillet`)
# ### Straight-line filleting
# For straight lines, the fillet solve is trivial and should converge in one iteration (as long as the fillet radius is shorter than both lines).

p_corner = SVector(0.0, 0.0)

lineA = LineSeg(SVector(1.0, -1.0), p_corner)
lineB = LineSeg(p_corner, SVector(1.0, 1.0))

plot(lineA, aspect_ratio = :equal)
plot!(lineB)

# The below gets summarized into `solve_fillet_arc()`

R = 0.25
trims, center = solve_fillet(R, (lineA, lineB))

trim_pt_A = lineA(trims[1])
trim_pt_B = lineB(trims[2])

vec_A = trim_pt_A - center
vec_B = trim_pt_B - center

θ₁ = mod(atan(vec_A[2], vec_A[1]), 2π)
θ₂ = mod(atan(vec_B[2], vec_B[1]), 2π)

fillet_arc = CircArc(center, R, θ₁, θ₂)

plot!(fillet_arc, linecolor = :blue)

# ### Fillet between two arcs

arcA = CircArc(SVector(-1.0, 0.0), 1.0, -π / 2, 0.0)
arcB = CircArc(SVector(0.0, 1.0), 1.0, -π / 2, 0.0)

plot(arcA, aspect_ratio = :equal)
plot!(arcB)

R = 0.25
trims, fillet_arc = solve_fillet_arc(R, (arcA, arcB))
plot!(fillet_arc, linecolor = :blue)

# ### Filleting a cusp

arcA = CircArc(SVector(-1.0, 0.0), 1.0, -π / 2, 0.0)
arcB = CircArc(SVector(1.0, 0.0), 1.0, 1.0 * π, 3 * π / 2)

plot(arcA, aspect_ratio = :equal)
plot!(arcB)

R = 0.2
trims, fillet_arc = solve_fillet_arc(R, (arcA, arcB))
plot!(fillet_arc, linecolor = :blue)


# ### Weird one

arcA = CircArc(SVector(-1.0, 0.0), 1.0, -π / 6, 0.0)
arcB = CircArc(SVector(-sqrt(3) / 2, -0.5), 1.0, π / 6, 7π / 6)

plot(arcA, aspect_ratio = :equal)
plot!(arcB)

R = 0.8
trims, fillet_arc = solve_fillet_arc(R, (arcA, arcB))
plot!(fillet_arc, linecolor = :blue)
