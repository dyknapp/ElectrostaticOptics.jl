# # Building a mesh

using ElectrostaticOptics
using StaticArrays
using Plots

# ### two circles

circle_mesh_A = circle(SVector(0.0, 2.0), 1.0)
circle_mesh_B = circle(SVector(5.0, 2.0), 1.0)

combined_mesh = concat_meshes(circle_mesh_A, circle_mesh_B)

plot(combined_mesh)

# ### Rectangles

rect = rectangle(
    SVector(10.0, 10.0),
    3.0, 4.0, π / 4.0
)

circles_and_rect = concat_meshes(combined_mesh, rect)

plot(circles_and_rect)

# ### Graph stuff
# Graphs.jl gives us some nice ways to ask questions about the mesh.
# How many connected components are there?  (There are two circles + 1 rectangle, so there should be 3).
using Graphs

count_connected_components(circles_and_rect.topology)


# ### Cycles
# How many cycles does the full graph have?  (Should have 3 again, circles and rectangle)
# The two two-element cycles are the circles, the rectangle is a four-element cycle.

simplecycles(circles_and_rect.topology)

# ### Degree of each vertex

for vert in vertices(circles_and_rect.topology)
    @assert ((indegree(circles_and_rect.topology, vert) == 1) && (outdegree(circles_and_rect.topology, vert) == 1))
end

# ### DOF assignment

assign_dofs()
