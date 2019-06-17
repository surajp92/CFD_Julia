clearconsole()

using CPUTime
using Printf
using ASTInterpreter2
using Plots

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

function restriction(nxf, nyf, nxc, nyc, r, sc)

    for j = 2:nyc for i = 2:nxc
        # grid index for fine grid for the same coarse point
        center = 4.0*r[2*i-1, 2*j-1]
        # E, W, N, S with respect to coarse grid point in fine grid
        grid = 2.0*(r[2*i-1, 2*j-1+1] + r[2*i-1, 2*j-1-1] +
                    r[2*i-1+1, 2*j-1] + r[2*i-1-1, 2*j-1])
        # NE, NW, SE, SW with respect to coarse grid point in fine grid
        corner = 1.0*(r[2*i-1+1, 2*j-1+1] + r[2*i-1+1, 2*j-1-1] +
                      r[2*i-1-1, 2*j-1+1] + r[2*i-1-1, 2*j-1-1])
        # restriction using trapezoidal rule
        sc[i,j] = (center + grid + corner)/16.0
    end end

    # restriction for boundary points bottom and top
    for j = 1:nxc+1
        # bottom boundary i = 1
        sc[1,j] = r[1, 2*j-1]
        # top boundary i = ny_coarse+1
        sc[nyc+1] = r[nyf+1, 2*j-1]
    end

    # restriction for boundary poinys left and right
    for i = 1:nyc+1
        # left boundary j = 1
        sc[i,1] = r[2*i-1,1]
        # right boundary nx_coarse+1
        sc[i,nxc+1] = r[2*i-1, nyf+1]
    end
end

function prolongation(nxc, nyc, nxf, nyf, unc, prol_fine)
    for j = 1:nyc for i = 1:nxc
        # direct injection at center point
        prol_fine[2*i-1, 2*j-1] = unc[i,j]
        # east neighnour on fine grid corresponding to coarse grid point
        prol_fine[2*i-1, 2*j-1+1] = 0.5*(unc[i,j] + unc[i,j+1])
        # north neighbout on fine grid corresponding to coarse grid point
        prol_fine[2*i-1+1, 2*j-1] = 0.5*(unc[i,j] + unc[i+1,j])
        # NE neighbour on fine grid corresponding to coarse grid point
        prol_fine[2*i-1+1, 2*j-1+1] = 0.25*(unc[i,j] + unc[i,j+1] +
                                            unc[i+1,j] + unc[i+1,j+1])

    end end

    # update boundary points
    for i = 1:nyc+1
        # left boundary j = 1
        prol_fine[2*i-1,1] = unc[i,1]
        # right boundary j = nx_fine+1
        prol_fine[2*i-1, nyf+1] = unc[i,nxc+1]
    end

    for j = 1:nxc+1
        #bottom boundary i = 1
        prol_fine[1,2*j-1] = unc[1,j]
        # top boundary i =  ny_fine+1
        prol_fine[nyf+1,2*j-1] = unc[nyc+1,j]
    end
end

function gauss_seidel_mg(nx, ny, dx, dy, f, un, V)

    rt = zeros(Float64, nx+1, ny+1)
    den = -2.0/dx^2 - 2.0/dy^2

    for iteration_count = 1:V
        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for j = 2:nx for i = 2:ny
            rt[i,j] = f[i,j] -
                     (un[i+1,j] - 2.0*un[i,j] + un[i-1,j])/dx^2 -
                     (un[i,j+1] - 2.0*un[i,j] + un[i,j-1])/dy^2

            un[i,j] = un[i,j] + rt[i,j]/den
        end end
    end
end

function mg(dx, dy, nx, ny, r, f, u_n, rms, v1, v2, v3, init_rms, max_iter, output)

    # create text file for writing residual history
    residual_plot = open("residual.csv", "w")
    write(residual_plot, "k"," ","rms"," ","rms/rms0"," \n")

    count = 0.0
    # compute initial residual
    compute_residual(nx, ny, dx, dy, f, u_n, r)
    # compute initial L-2 norm
    rms = compute_l2norm(nx, ny, r)

    init_rms = rms
    println("0", " ", rms, " ", rms/init_rms)


    #allocate memory for grid size at different levels
    lnx = zeros(Int64, 2)
    lny = zeros(Int64, 2)
    ldx = zeros(Float64, 2)
    ldy = zeros(Float64, 2)
    lnx[1] = nx
    lny[1] = ny
    lnx[2] = Int64(lnx[1]/2)
    lny[2] = Int64(lny[1]/2)
    ldx[1] = dx
    ldy[1] = dy
    ldx[2] = ldx[1]*2
    ldy[2] = ldy[1]*2

    # allocate matrix for storage at fine level
    # residual at fine level is already defined at global level
    prol_fine = zeros(Float64, lnx[1]+1, lny[1]+1)

    # allocate matrix for storage at coarse levels
    fc = zeros(Float64, lnx[2]+1, lny[2]+1)
    unc = zeros(Float64, lnx[2]+1, lny[2]+1)

    # start main iteration loop
    for iteration_count = 1:max_iter
        # call relaxation on fine grid and compute the numerical solution
        # for fixed number of iteration

        gauss_seidel_mg(lnx[1], lny[1], dx, dy, f, u_n, v1)

        # check for convergence only for finest grid
        # compute the residual and L2 norm

        compute_residual(nx, ny, dx, dy, f, u_n, r)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, r)
        # write results only for finest residual
        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/init_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/init_rms)

        if (rms/init_rms) <= tolerance
                break
        end

        # restrict the residual from fine level to coarse level

        restriction(lnx[1], lny[1], lnx[2], lny[2], r, fc)

        # set solution zero on coarse grid
        unc[:,:] = zeros(lnx[2]+1, lny[2]+1)

        # solve on the coarsest level and relax V3 times
        gauss_seidel_mg(lnx[2], lny[2], ldx[2], ldy[2],
                        fc, unc, v3)

        # prolongate solution from coarse level to fine level
        prolongation(lnx[2], lny[2], lnx[1], lny[1], unc, prol_fine)

        # correct the solution on fine level
        for j = 2:lnx[1] for i = 2:lny[1]
                u_n[i,j] = u_n[i,j] + prol_fine[i,j]
        end end

        # relax v2 times
        gauss_seidel_mg(lnx[1], lny[1], dx, dy, f, u_n, v2)
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

write(field_initial, "x y f un ue \n")
write(field_final, "x y f un ue e \n")

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

v1 = 2 # relaxation
v2 = 2 # prolongation
v3 = 2 # coarsest level
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
    # u_e[i,j] = sin(3.0*x[i]) + cos(2.0*y[j])
    # f[i,j] = -9.0*sin(3.0*x[i]) -4.0*cos(2.0*y[j])
    # u_n[i,j] = 0.0
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

mg(dx, dy, nx, ny, r, f, u_n, rms, v1, v2, v3, init_rms, max_iter, output)

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
close(output)
