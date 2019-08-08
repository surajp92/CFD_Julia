#clearconsole()

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
    rms = 0.0
    # println(residual)
    for j = 1:ny+1 for i = 1:nx+1
        rms = rms + r[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx+1)*(ny+1)))
    return rms
end

#-----------------------------------------------------------------------------#
# Fast poisson solver for periodic domain
#-----------------------------------------------------------------------------#
function fps(nx,ny,dx,dy,f)
    eps = 1.0e-6

    kx = Array{Float64}(undef,nx)
    ky = Array{Float64}(undef,ny)

    data = Array{Complex{Float64}}(undef,nx,ny)
    data1 = Array{Complex{Float64}}(undef,nx,ny)
    e = Array{Complex{Float64}}(undef,nx,ny)

    u = Array{Complex{Float64}}(undef,nx,ny)

    aa = -2.0/(dx*dx) - 2.0/(dy*dy)
    bb = 2.0/(dx*dx)
    cc = 2.0/(dy*dy)

    #wave number indexing
    hx = 2.0*pi/nx

    for i = 1:Int64(nx/2)
        kx[i] = hx*(i-1.0)
        kx[i+Int64(nx/2)] = hx*(i-Int64(nx/2)-1)
    end
    kx[1] = eps

    ky = kx

    for i = 1:nx
        for j = 1:ny
            data[i,j] = complex(f[i,j],0.0)
        end
    end

    e = fft(data)
    e[1,1] = 0.0
    for i = 1:nx
        for j = 1:ny
            data1[i,j] = e[i,j]/(aa + bb*cos(kx[i]) + cc*cos(ky[j]))
        end
    end

    u = real(ifft(data1))

    return u
end

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
function numerical(nx,ny,nt,dx,dy,dt,re,wn,ns)

    wt = Array{Float64}(undef, nx+2, ny+2) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx+2, ny+2)
    ut = Array{Float64}(undef, nx+1, ny+1)

    m = 1 # record index
    freq = Int64(nt/ns)

    for k = 1:nt
        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,r)

        for i = 2:nx+1 for j = 2:ny+1
            wt[i,j] = wn[i,j] + dt*r[i,j]
        end end

        # periodic BC
        wt[nx+2,:] = wt[2,:]
        wt[:,ny+2] = wt[:,2]

        # ghost points
        wt[1,:] = wt[nx+1,:]
        wt[:,1] = wt[:,ny+1]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,r)

        for i = 2:nx+1 for j = 2:ny+1
            wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]
        end end

        # periodic BC
        wt[nx+2,:] = wt[2,:]
        wt[:,ny+2] = wt[:,2]

        # ghost points
        wt[1,:] = wt[nx+1,:]
        wt[:,1] = wt[:,ny+1]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,r)

        for i = 2:nx+1 for j = 2:ny+1
            wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]
        end end

        # periodic BC
        wn[nx+2,:] = wn[2,:]
        wn[:,ny+2] = wn[:,2]

        # ghost points
        wn[1,:] = wn[nx+1,:]
        wn[:,1] = wn[:,ny+1]

        if (mod(k,freq) == 0)
            println(k)
            ut = wn[2:nx+2,2:ny+2]
            field_final = open(string("vm",string(m),".txt"), "w");
            for j = 1:ny+1 for i = 1:nx+1
                write(field_final, string(x[i]), " ",string(y[j]), " ", string(ut[i,j]), " \n")
            end end
            m = m+1
            close(field_final)
        end
    end

    return wn[2:nx+2,2:ny+2]
end

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -J(w,ψ) + ν ∇^2(w)
#-----------------------------------------------------------------------------#
function rhs(nx,ny,dx,dy,re,w,r)

    # compute streamfunction from vorticity
    s = Array{Float64}(undef, nx+2, ny+2)
    f = Array{Float64}(undef, nx, ny)
    f = -w[2:nx+1,2:ny+1]

    s[2:nx+1,2:ny+1] = fps(nx,ny,dx,dy,f)

    # periodic BC
    s[nx+2,:] = s[2,:]
    s[:,ny+2] = s[:,2]

    # ghost points
    s[1,:] = s[nx+1,:]
    s[:,1] = s[:,ny+1]

    # Arakawa numerical scheme for Jacobian
    aa = 1.0/(re*dx*dx)
    bb = 1.0/(re*dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0

    for i = 2:nx+1 for j = 2:ny+1
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

# initial condition for vortex merger problem
function vm_ic(nx,ny,x,y,w)
    sigma = pi
    xc1 = pi-pi/4.0
    yc1 = pi
    xc2 = pi+pi/4.0
    yc2 = pi

    for i = 2:nx+2 for j = 2:ny+2
        w[i,j] = exp(-sigma*((x[i-1]-xc1)^2 + (y[j-1]-yc1)^2)) +
                 exp(-sigma*((x[i-1]-xc2)^2 + (y[j-1]-yc2)^2))
    end end
end


#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
nx = 128
ny = 128

x_l = 0.0
x_r = 2.0*pi
y_b = 0.0
y_t = 2.0*pi

dx = (x_r-x_l)/nx
dy = (y_t-y_b)/ny

dt = 0.01
tf = 20.0
nt = tf/dt
re = 1000.0
ns = 10

x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)

for i = 1:nx+1
    x[i] = dx*(i-1)
end
for i = 1:ny+1
    y[i] = dy*(i-1)
end

wn = Array{Float64}(undef, nx+2, ny+2)
un = Array{Float64}(undef, nx+1, ny+1)
un0 = Array{Float64}(undef, nx+1, ny+1)
ue = Array{Float64}(undef, nx+1, ny+1)
uerror = Array{Float64}(undef, nx+1, ny+1)

time = 0.0

vm_ic(nx,ny,x,y,wn)
# ghost points
wn[1,:] = wn[nx+1,:]
wn[:,1] = wn[:,ny+1]

wn[nx+2,:] = wn[2,:]
wn[:,ny+2] = wn[:,2]

un0 = wn[2:nx+2,2:ny+2]

field_final = open("vm0.txt", "w");
for j = 1:ny+1 for i = 1:nx+1
    write(field_final, string(x[i]), " ",string(y[j]), " ", string(un0[i,j]), " \n")
end end

val, t, bytes, gctime, memallocs = @timed begin
un = numerical(nx,ny,nt,dx,dy,dt,re,wn,ns)
end

print("CPU Time = ", t);

time = tf

field_final = open("field_final.txt", "w");
for j = 1:ny+1 for i = 1:nx+1
    write(field_final, string(x[i]), " ",string(y[j]), " ", string(un[i,j]), " \n")
end end
close(field_final)
