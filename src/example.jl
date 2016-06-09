using Gadfly
include("solver.jl")

c = solver()

# Plot initial concentration
plot_concentration(c, 0)
# Plot final concentration
plot_concentration(c, 36500)
