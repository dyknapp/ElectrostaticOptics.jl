# # Lagrange Basis

using ElectrostaticOptics
using StaticArrays
using Plots

# ## Different orders
# ### 0th order

lb0 = LagrangeBasis(Float64, 0)
plot(lb0)

# ### 1st order

lb1 = LagrangeBasis(Float64, 1)
plot(lb1)

# ### 2nd order

lb2 = LagrangeBasis(Float64, 2)
plot(lb2)

# ### 3rd order

lb3 = LagrangeBasis(Float64, 3)
plot(lb3)

# ### 4th order

lb4 = LagrangeBasis(Float64, 4)
plot(lb4)
