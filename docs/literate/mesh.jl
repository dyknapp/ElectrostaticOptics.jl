# # Building a mesh

using ElectrostaticOptics
using StaticArrays
using Plots

# ### two circles

circle_mesh_A = circle(SVector(0.0, 2.0), 1.0)
circle_mesh_B = circle(SVector(5.0, 2.0), 1.0)

combined_mesh = concat_meshes(circle_mesh_A, circle_mesh_B)

plot(combined_mesh)


# ### Graph stuff
# Graphs.jl gives us some nice ways to ask questions about the mesh.
# How many connected components are there?  (There are two circles, so there should be 2).
using Graphs

Graphs.count_connected_components(combined_mesh.topology)


# ### Rectangles

rect = rectangle(
    SVector(10.0, 10.0),
    3.0, 4.0, π / 4.0
)

circles_and_rect = concat_meshes(combined_mesh, rect)

plot(circles_and_rect)
