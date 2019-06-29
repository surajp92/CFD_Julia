using CPUTime
using Printf
using Plots
#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
function compute_l2norm(nx,r)
    rms = 0.0
    for i = 2:nx
        rms = rms + r[i]^2
    end
    rms = sqrt(rms/((nx-1)))
    return rms
end

#-----------------------------------------------------------------------------#
# Solution to tridigonal system using Thomas algorithm
#-----------------------------------------------------------------------------#
function tdma(a,b,c,r,x,s,e)
    for i = s+1:e
        b[i] = b[i] - a[i]*(c[i-1]/b[i-1])
        r[i] = r[i] - a[i]*(r[i-1]/b[i-1])
    end

    x[e] = r[e]/b[e]

    for i = e-1:-1:s
        x[i] = (r[i] - c[i]*x[i+1])/b[i]
    end
end

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 5th-order Compact WENO scheme for spatial terms
#-----------------------------------------------------------------------------#
function numerical(nx,nt,dx,dt,x,u,α)
    un = Array{Float64}(undef, nx+1) # numerical solsution at every time step
    ut = Array{Float64}(undef, nx+1) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx)

    k = 1 # record index

    for i = 1:nx+1
        un[i] = -sin(pi*x[i])
        u[i,k] = un[i] # store solution at t=0
    end

    # dirichlet boundary condition
    un[1] = 0.0
    un[nx+1] = 0.0

    # dirichlet boundary condition for temporary array
    ut[1] = 0.0
    ut[nx+1] = 0.0

    for j = 2:nt+1
        rhs(nx,dx,un,r,α)

        for i = 2:nx
            ut[i] = un[i] + dt*r[i]
        end

        rhs(nx,dx,ut,r,α)

        for i = 2:nx
            ut[i] = 0.75*un[i] + 0.25*ut[i] + 0.25*dt*r[i]
        end

        rhs(nx,dx,ut,r,α)

        for i = 2:nx
            un[i] = (1.0/3.0)*un[i] + (2.0/3.0)*ut[i] + (2.0/3.0)*dt*r[i]
        end

        k = k+1
        u[:,k] = un[:]
    end
end

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -u∂u/∂x
#-----------------------------------------------------------------------------#
function rhs(nx,dx,u,r,α)
    for i = 2:nx
        r[i] = α*(u[i+1] - 2.0*u[i] + u[i-1])/(dx*dx)
    end
end

#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
x_l = -1.0
x_r = 1.0
dx = 0.025
nx = Int64((x_r-x_l)/dx)

dt = 0.0025
t = 1.0
nt = Int64(t/dt)

α = 1/(pi*pi)

u = Array{Float64}(undef, nx+1, nt+1)
x = Array{Float64}(undef, nx+1)
u_e = Array{Float64}(undef, nx+1)
error = Array{Float64}(undef, nx+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)  # location of each grid point
    u_e[i] = -exp(-t)*sin(pi*x[i]) # initial condition @ t=0
end

numerical(nx,nt,dx,dt,x,u,α)

u_error = zeros(nx+1)
u_error = u[:,nt+1] - u_e
rms_error = compute_l2norm(nx, u_error)
max_error = maximum(abs.(u_error))

# create output file for L2-norm
output = open("output.txt", "w");
write(output, "Error details: \n");
write(output, "L-2 Norm = ", string(rms_error), " \n");
write(output, "Maximum Norm = ", string(max_error), " \n");

# create text file for final field
field_final = open("field_final.csv", "w");
write(field_final, "x"," ", "ue", " ", "un", " ", "uerror" ," \n")

for i = 1:nx+1
    write(field_final, @sprintf("%.16f",x[i])," ",@sprintf("%.16f", u_e[i])," ",
          @sprintf("%.16f", u[i,nt+1])," ",@sprintf("%.16f", u_error[i])," \n")
end

close(field_final)
close(output);
