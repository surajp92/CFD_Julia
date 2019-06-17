clearconsole()

using CPUTime
using Printf
using Plots
using FFTW

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

nx = 128
ny = 128

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

dx = (x_r-x_l)/nx
dy = (y_t-y_b)/ny

x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)
ue = Array{Float64}(undef,nx+1,ny+1)
f = Array{Float64}(undef,nx+1,ny+1)
un = Array{Float64}(undef,nx+1,ny+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)
end
for i = 1:ny+1
    y[i] = y_b + dy*(i-1)
end

# given exact solution
km = 16
c1 = (1.0/km)^2
c2 = -8.0*pi*pi

for j = 1:ny+1
    for i = 1:nx+1
        # ue[i,j] = sin(3.0*2.0*pi*x[i]) + cos(2.0*2.0*pi*y[j])
        # f[i,j] = -9.0*4.0*pi^2*sin(3.0*2.0*pi*x[i]) - 4.0*4.0*pi^2*cos(2.0*2.0*pi*y[j])

        # ue[i,j] = cos(2.0*pi*x[i]) + cos(2.0*pi*y[j])
        # f[i,j] = -4.0*pi*pi*ue[i,j]

        ue[i,j] = sin(2.0*pi*x[i])*sin(2.0*pi*y[j]) +
               c1*sin(km*2.0*pi*x[i])*sin(km*2.0*pi*y[j])

        f[i,j] = c2*sin(2.0*pi*x[i])*sin(2.0*pi*y[j]) +
                 c2*sin(km*2.0*pi*x[i])*sin(km*2.0*pi*y[j])

        un[i,j] = 0.0
    end
end

val, t, bytes, gctime, memallocs = @timed begin

un[1:nx,1:ny] = fps(nx,ny,dx,dy,f)

end

# Periodic boundary condition
un[nx+1,:] = un[1,:]
un[:,ny+1] = un[:,1]

uerror = zeros(nx+1, ny+1)
rms_error = 0.0

uerror = un - ue

rms_error = compute_l2norm(nx, ny, uerror)
max_error = maximum(abs.(uerror))

println("Error details:");
println("L-2 Norm = ", rms_error);
println("Maximum Norm = ", max_error);
print("CPU Time = ", t);

p1 = contour(x, y, transpose(ue), fill=true)
p2 = contour(x, y, transpose(un), fill=true)
p3 = plot(p1,p2, size = (1000, 400))
savefig(p3,"contourd.pdf")
