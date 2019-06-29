clearconsole()

using CPUTime
using Printf
using Plots
using ASTInterpreter2

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

x_l = -1.0
x_r = 1.0
dx = 0.025
nx = Int64((x_r-x_l)/dx)

dt = 0.0025
t = 1.0
nt = Int64(t/dt)

α = 1/(pi*pi)

x = Array{Float64}(undef, nx+1)
u_e = Array{Float64}(undef, nx+1)
un = Array{Float64}(undef, nt+1, nx+1)
error = Array{Float64}(undef, nx+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)  # location of each grid point
    un[1,i] = -sin(pi*x[i]) # initial condition @ t=0
    u_e[i] = -exp(-t)*sin(pi*x[i]) # initial condition @ t=0
end

un[1,1] = 0.0
un[1,nx+1] = 0.0

beta = α*dt/(dx*dx)

for k = 2:nt+1
    for i = 2:nx
        un[k,i] = un[k-1,i] + beta*(un[k-1,i+1] -
                                2.0*un[k-1,i] + un[k-1,i-1])
    end
    un[k,1] = 0.0 # boundary condition at x = -1
    un[k,nx+1] = 0.0 # boundary condition at x = -1
end

# compute L2 norm of the error
u_error = un[nt+1,:] - u_e
rms_error = compute_l2norm(nx,u_error)
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
          @sprintf("%.16f", un[nt+1,i])," ",@sprintf("%.16f", u_error[i])," \n")
end

close(field_final)
close(output);
