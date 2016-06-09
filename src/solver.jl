"""
The code we went over. The formulas are derived in the notes I sent in the
email. If things start blowing up again, try making the time step (k) smaller.
This would be cleaner if you make a function for the main solver and use as the
space step and time step as inputs, i.e.

function solver(h,k)
    ...
    return c
end

That can be the first thing to do as you learn Julia :)

Also, remember we played around with the dispersion gradient. Make sure that
you are okay with the formulation before you trust the output answers.
"""


################################################################################
# Model Functions
function velocity(x, t)
    h0 = 1.0
    k = 2.0
    beta = 0.122
    t0 = 0.52
    ne = 0.2
    ir = -0.01
    vel = h0*k/ne*beta*exp(-1.0*beta*x)*sqrt(2.0)*sin(2.0*pi*t/t0 - x*beta + pi/4.0) + k/ne*ir
    return vel
end

function dispersion(x, t)
    alpha = 10.0
    disp = alpha*abs(velocity(x, t))
    return disp
end

function velocity_grad(x, t)
    h0 = 1.0
    k = 2.0
    beta = 0.122
    t0 = 0.52
    ne = 0.2
    grad = -1.0*sqrt(2.0)*beta^2*h0*k/ne*exp(-1.0*beta*x)*sin(-beta*x + 2*pi*t/t0 + pi/4) -
    sqrt(2.0)*beta^2*h*k/ne*exp(-1.0*beta*x)*cos(-1.0*beta*x + 2*pi*t/t0 + pi/4)
    return grad
end

#function dispersion_grad(x, t)
#    alpha = 10.0
#    disp_grad = alpha*velocity(x,t)/abs(velocity(x,t))*velocity_grad(x,t)
#    return disp_grad
#end

function dispersion_grad(x, t)
    alpha = 10.0
    disp_grad = alpha*velocity_grad(x,t)
    return disp_grad
end


function initial_condition(x)
    if x > 160
        ic = 0.0
    else
        ic = -1.0/80.0*abs(x - 80.0) + 1.0
    end
    return ic
end


################################################################################
# Main finite difference solver

vmax = 1.2 # ft/day
h = 1.0 # ft, grid step size
#k = h/vmax # days, time step size
k = 0.01
#L = 10.0 # ft, domain
#T = 10.0 # days, end time
#n = ceil(Int, L/h + 1) # number of space grid points
#m = ceil(Int, T/k + 1) # number of time grid points
n = 200
m = 36500
c = Array{Float64}(m,n)

# initial conditions
for i = 1:n
    c[1,i] = initial_condition((i-1)*h)
end

# boundary conditions
#for j = 1:m
#    if velocity(0, m-1) > 0
#        # Do condition 1
#    else
#        # Do condition 2
#    end
#end

for j = 1:m
    c[j,1] = 0.0
    c[j,n] = 0.0
end

# Main solver
for i = 2:m
    for j = 2:n-1
        x = (j-1)*h
        t = (i-1)*k
        coef1 = k*dispersion(x,t)/h^2 - k/(2*h)*(dispersion_grad(x,t) - velocity(x,t))
        coef2 = 1 - k*velocity_grad(x,t) - 2*k*dispersion(x,t)/h^2
        coef3 = k*dispersion(x,t)/h^2 + k/(2*h)*(dispersion_grad(x,t) - velocity(x,t))
        c[i,j] = coef1*c[i-1,j-1] + coef2*c[i-1,j] + coef3*c[i-1,j+1]
    end
end
################################################################################
# Plotting code: This will plot the initial condition and the concentration after
# 1 year. Loading Gadfly is slow the first time you use it and the first plot
# is slow the first time you use it.
using Gadfly # This takes a while to load the first time
x_plot = linspace(0,200,n) # n uniformly spaced plot points from 0 ft to 200 ft
plot(x=x_plot, y=c[1,:], Geom.point) # Initial condition
plot(x=x_plot, y=c[m,:], Geom.point) # Final time
