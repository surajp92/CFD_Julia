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

function conjugate_gradient(dx, dy, nx, ny, r, f, u_n, rms,
                      init_rms, max_iter, tolerance, tiny, output)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    #write(residual_plot, "k"," ","rms"," ","rms/rms0"," \n")

    count = 0.0

    compute_residual(nx, ny, dx, dy, f, u_n, r)

    rms = compute_l2norm(nx, ny, r)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", initial_rms, " ", rms/initial_rms)
    # allocate the matric for direction and set the initial direction (conjugate vector)
    p = zeros(Float64, nx+1, ny+1)

    # asssign conjugate vector to initial residual
    for j = 1:ny+1 for i = 1:nx+1
        p[i,j] = r[i,j]
    end end

    del_p    = zeros(Float64, nx+1, ny+1)

    # start calculation
    for iteration_count = 1:max_iter

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx
            del_p[i,j] = (p[i+1,j] - 2.0*p[i,j] + p[i-1,j])/(dx^2) +
                         (p[i,j+1] - 2.0*p[i,j] + p[i,j-1])/(dy^2)
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + r[i,j]*r[i,j]
            bb = bb + del_p[i,j]*p[i,j]
        end end
        # cc = <r,r>/<d,p>
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of conjugate vector
        for j = 2:ny for i = 2:nx
            u_n[i,j] = u_n[i,j] + cc*p[i,j]
        end end

        # bb = <r,r> = aa (calculated in previous loop)
        bb = aa
        aa = 0.0

        # update the residual by removing some component of previous residual
        for j = 2:ny for i = 2:nx
            r[i,j] = r[i,j] - cc*del_p[i,j]
            aa = aa + r[i,j]*r[i,j]
        end end
        # cc = <r-cd, r-cd>/<r,r>
        cc = aa/(bb+tiny)

        # update the conjugate vector
        for j = 1:ny for i = 1:nx
            p[i,j] = r[i,j] + cc*p[i,j]
        end end

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, r)

        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
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
tiny = Float64(1.0e-16)

# create output file for L2-norm
output = open("output.txt", "w");
write(output, "Residual details: \n");
# create text file for initial and final field
field_initial = open("field_initial.txt", "w")
field_final = open("field_final.txt", "w")

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

    u_e[i,j] = sin(2.0*pi*x[i]) * sin(2.0*pi*y[j]) +
               c1*sin(16.0*pi*x[i]) * sin(16.0*pi*y[j])

    f[i,j] = 4.0*c2*sin(2.0*pi*x[i]) * sin(2.0*pi*y[j]) +
                  c2*sin(16.0*pi*x[i]) * sin(16.0*pi*y[j])

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
    write(field_initial, string(x[i]), " ",string(y[j]), " ", string(f[i,j]),
          " ", string(u_n[i,j]), " ", string(u_e[i,j]), " \n")
end end
val, t, bytes, gctime, memallocs = @timed begin

conjugate_gradient(dx, dy, nx, ny, r, f, u_n, rms,
                      init_rms, max_iter, tolerance, tiny, output)

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
        write(field_final, string(x[i]), " ",string(y[j]), " ", string(f[i,j]),
              " ", string(u_n[i,j]), " ", string(u_e[i,j]), " \n")
end end

close(field_initial)
close(field_final)
close(output);
