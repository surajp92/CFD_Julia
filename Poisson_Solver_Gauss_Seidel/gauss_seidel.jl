clearconsole()

using CPUTime
using Printf
using ASTInterpreter2

function compute_residual(nx, ny, dx, dy, f, u_n, r)

    for j = 2:ny for i = 2:nx
        d2udx2 = (u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j])/(dx^2)
        d2udy2 = (u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])/(dy^2)
        r[i,j] = f[i,j]  - d2udx2 - d2udy2
    end end

end

function compute_l2norm(nx, ny, r)

    rms = 0.0
    # println(residual)
    for j = 2:ny for i = 2:nx
        rms = rms + r[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx-1)*(ny-1)))
    return rms
end

function gauss_seidel(dx, dy, nx, ny, r, f, u_n, rms, init_rms, max_iter,
                      tolerance, output)
    # create text file for writing residual history
    residual_plot = open("residual.csv", "w")
    write(residual_plot, "k"," ","rms"," ","rms/rms0"," \n")

    count = 0.0

    compute_residual(nx, ny, dx, dy, f, u_n, r)

    rms = compute_l2norm(nx, ny, r)

    init_rms = rms
    iter_count = 0
    println(iter_count, " ", init_rms, " ", rms/init_rms)

    den = -2.0/dx^2 - 2.0/dy^2
    for iter_count = 1:max_iter

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        # residual = f + λ^2u - ∇^2u
        for j = 2:ny for i = 2:nx
            d2udx2 = (u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j])/(dx^2)
            d2udy2 = (u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])/(dy^2)
            r[i,j] = f[i,j] - d2udx2 - d2udy2

            u_n[i,j] = u_n[i,j] + r[i,j]/den
        end end

        compute_residual(nx, ny, dx, dy, f, u_n, r)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, r)

        write(residual_plot, string(iter_count), " ",string(rms), " ", string(rms/init_rms)," \n");
        count = iter_count

        println(iter_count, " ", rms, " ", rms/init_rms)

        if (rms/init_rms) <= tolerance
            break
        end
    end

    write(output, "L-2 Norm = ", string(rms), " \n");
    write(output, "Maximum Norm = ", string(maximum(abs.(r))), " \n");
    write(output, "Iterations = ", string(count), " \n");
    close(residual_plot)
end

nx = Int64(128)
ny = Int64(128)
tolerance = Float64(1.0e-12)
max_iter = Int64(100000)

# create output file for L2-norm
output = open("output.txt", "w");
write(output, "Residual details: \n");

# create text file for initial and final field
field_initial = open("field_initial.csv", "w");
field_final = open("field_final.csv", "w");

#write(field_initial, "x"," ","y"," ","f"," ","un"," ","ue"," ","e", "\n")
#write(field_final, "x"," ","y"," ","f"," ","un"," ","ue"," ","e", "\n")

write(field_initial, "x y f un ue \n")
write(field_final, "x y f un ue e \n")

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

dx = (x_r - x_l)/nx
dy = (y_t - y_b)/ny

# allocate array for x and y position of grids, exact solution and source term
x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)
u_e = Array{Float64}(undef, nx+1, ny+1)
f = Array{Float64}(undef, nx+1, ny+1)
u_n = Array{Float64}(undef, nx+1, ny+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)
end
for i = 1:ny+1
    y[i] = y_b + dy*(i-1)
end

c1 = (1.0/16.0)^2
c2 = -2.0*pi*pi

for i = 1:nx+1 for j = 1:ny+1
    f[i,j] = c2 * sin(pi*x[i]) * sin(pi*y[j]) +
                  c2*sin(16.0*pi*x[i]) * sin(16.0*pi*y[j])

    u_e[i,j] = sin(pi*x[i]) * sin(pi*y[j]) +
               c1*sin(16.0*pi*x[i]) * sin(16.0*pi*y[j])

    u_n[i,j] = 0.0
end end

u_n[:,1] = u_e[:,1]
u_n[:, ny+1] = u_e[:, ny+1]

u_n[1,:] = u_e[1,:]
u_n[nx+1,:] = u_e[nx+1,:]

r = zeros(Float64, nx+1, ny+1)
init_rms = 0.0
rms = 0.0

for j = 1:ny+1 for i = 1:nx+1
    write(field_initial, @sprintf("%.16f",x[i])," ", @sprintf("%.16f", y[j]), " ",
          @sprintf("%.16f", f[i,j])," ", @sprintf("%.16f", u_n[i,j])," ",
          @sprintf("%.16f", u_e[i,j]), " \n")
end end
val, t, bytes, gctime, memallocs = @timed begin

gauss_seidel(dx, dy, nx, ny, r, f, u_n, rms, init_rms, max_iter, tolerance, output)

end
u_error = zeros(nx+1, ny+1)
rms_error = 0.0

u_error = u_n - u_e
rms_error = compute_l2norm(nx, ny, u_error)
max_error = maximum(abs.(u_error))

println("Error details:");
println("L-2 Norm = ", rms_error);
println("Maximum Norm = ", max_error);
print("CPU Time = ", t);

write(output, "Error details: \n");
write(output, "L-2 Norm = ", string(rms_error), " \n");
write(output, "Maximum Norm = ", string(max_error), " \n");
write(output, "CPU Time = ", string(t), " \n");

for j = 1:ny+1 for i = 1:nx+1
    write(field_final, @sprintf("%.16f",x[i])," ", @sprintf("%.16f", y[j]), " ",
          @sprintf("%.16f", f[i,j])," ", @sprintf("%.16f", u_n[i,j])," ",
          @sprintf("%.16f", u_e[i,j])," ", @sprintf("%.16f",(u_error[i,j]))," \n")
end end

close(field_initial)
close(field_final)
close(output);
