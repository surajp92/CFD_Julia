clearconsole()

using CPUTime
using Printf
using Plots
using FFTW


font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
function compute_l2norm(nx, ny, r)
    norm = 0.0
    # println(residual)
    for j = 1:ny+1 for i = 1:nx+1
        norm = norm + r[i,j]^2
    end end
    # println(rms)
    norm = sqrt(norm/((nx+1)*(ny+1)))
    return norm
end

#-----------------------------------------------------------------------------#
# Fast poisson solver for homozeneous Dirichlet domain
#-----------------------------------------------------------------------------#
function fps_sine(nx,ny,dx,dy,f)

    data = Array{Complex{Float64}}(undef,nx-1,ny-1)
    data1 = Array{Complex{Float64}}(undef,nx-1,ny-1)
    e = Array{Complex{Float64}}(undef,nx-1,ny-1)

    u = Array{Complex{Float64}}(undef,nx-1,ny-1)

    for i = 1:nx-1
        for j = 1:ny-1
            data[i,j] = f[i+1,j+1]
        end
    end

    e = FFTW.r2r(data,FFTW.RODFT00)

    for i = 1:nx-1
        for j = 1:ny-1
            alpha = (2.0/(dx*dx))*(cos(pi*i/nx) - 1.0) +
                    (2.0/(dy*dy))*(cos(pi*j/ny) - 1.0)
            data1[i,j] = e[i,j]/alpha
        end
    end

    u = FFTW.r2r(data1,FFTW.RODFT00)/((2*nx)*(2*ny))

    return u
end

function bc(nx,ny,dx,dy,w,s)
    # first order approximation
    # boundary condition for vorticity (Hoffmann) left and right
    for j = 1:ny+1
        w[1,j] = -2.0*s[2,j]/(dx*dx)
        w[nx+1,j]= -2.0*s[nx,j]/(dx*dx)
    end

    # boundary condition for vorticity (Hoffmann) bottom and top
    for i = 1:nx+1
        w[i,1] = -2.0*s[i,2]/(dy*dy)
        w[i,ny+1]= -2.0*s[i,ny]/(dy*dy) - 2.0/dy
    end
end

function bc2(nx,ny,dx,dy,w,s)
    # second order approximation
    # boundary condition for vorticity (Jensen) left and right
    for j = 1:ny+1
        w[1,j] = (-4.0*s[2,j]+0.5*s[3,j])/(dx*dx)
        w[nx+1,j]= (-4.0*s[nx,j]+0.5*s[nx-1,j])/(dx*dx)
    end

    # boundary condition for vorticity (Jensen) bottom and top
    for i = 1:nx+1
        w[i,1] = (-4.0*s[i,2]+0.5*s[i,3])/(dy*dy)
        w[i,ny+1]= (-4.0*s[i,ny]+0.5*s[i,ny-1])/(dy*dy) - 3.0/dy
    end
end

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
function numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms)

    wt = Array{Float64}(undef, nx+1, ny+1) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx+1, ny+1) # right hand side
    sp = Array{Float64}(undef, nx+1, ny+1) # old streamfunction

    for k = 1:nt

        for i = 1:nx+1 for j = 1:ny+1
            sp[i,j] = sn[i,j]
        end end

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,sn,r)

        for i = 2:nx for j = 2:ny
            wt[i,j] = wn[i,j] + dt*r[i,j]
        end end
        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction from vorticity
        sn[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,-wt)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for i = 2:nx for j = 2:ny
            wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]
        end end
        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction from vorticity
        sn[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,-wt)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for i = 2:nx for j = 2:ny
            wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]
        end end
        bc2(nx,ny,dx,dy,wn,sn)

        # compute streamfunction from vorticity
        sn[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,-wn)

        rms[k] = 0.0
        for i = 1:nx+1 for j = 1:ny+1
            rms[k] = rms[k] + (sn[i,j] - sp[i,j])^2
        end end

        rms[k] = sqrt(rms[k]/((nx+1)*(ny+1)))
        println(k, " ", rms[k])

    end
end

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -J(w,ψ) + ν ∇^2(w)
#-----------------------------------------------------------------------------#
function rhs(nx,ny,dx,dy,re,w,s,r)
    # Arakawa numerical scheme for Jacobian
    aa = 1.0/(re*dx*dx)
    bb = 1.0/(re*dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0

    for i = 2:nx for j = 2:ny
        j1 = gg*((w[i+1,j]-w[i-1,j])*(s[i,j+1]-s[i,j-1]) -
                 (w[i,j+1]-w[i,j-1])*(s[i+1,j]-s[i-1,j]))

        j2 = gg*(w[i+1,j]*(s[i+1,j+1]-s[i+1,j-1]) -
                 w[i-1,j]*(s[i-1,j+1]-s[i-1,j-1]) -
                 w[i,j+1]*(s[i+1,j+1]-s[i-1,j+1]) +
                 w[i,j-1]*(s[i+1,j-1]-s[i-1,j-1]))

        j3 = gg*(w[i+1,j+1]*(s[i,j+1]-s[i+1,j]) -
                 w[i-1,j-1]*(s[i-1,j]-s[i,j-1]) -
            	 w[i-1,j+1]*(s[i,j+1]-s[i-1,j]) +
            	 w[i+1,j-1]*(s[i+1,j]-s[i,j-1]))

        jac = (j1+j2+j3)*hh

        #Central difference for Laplacian
        r[i,j] = -jac + aa*(w[i+1,j]-2.0*w[i,j]+w[i-1,j]) +
                        bb*(w[i,j+1]-2.0*w[i,j]+w[i,j-1])
        end end
end


#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
nx = 64
ny = 64

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

dx = (x_r-x_l)/nx
dy = (y_t-y_b)/ny

dt = 0.001
tf = 10.0
nt = Int64(tf/dt)
re = 100.0

x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)
rms = Array{Float64}(undef, nt)

for i = 1:nx+1
    x[i] = dx*(i-1)
end
for i = 1:ny+1
    y[i] = dy*(i-1)
end

wn = Array{Float64}(undef, nx+1, ny+1)
sn = Array{Float64}(undef, nx+1, ny+1)

time = 0.0

for i = 1:nx+1 for j = 1:ny+1
    wn[i,j] = 0.0 # initial condition
    sn[i,j] = 0.0 # initial streamfunction
end end

numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms)

field_final = open("field_final.txt", "w");
#write(field_final, "x y wn sn \n")

residual_plot = open("residual_plot.txt", "w");
#write(residual_plot, "n res \n")
t = Array(dt:dt:tf)

for n = 1:nt
    write(residual_plot, string(n), " ",string(rms[n])," \n");
end
close(residual_plot)

for j = 1:ny+1 for i = 1:nx+1
    write(field_final, string(x[i]), " ",string(y[j]), " ", string(wn[i,j]),
          " ", string(sn[i,j]), " \n")
end end
close(field_final)
