clearconsole()

using CPUTime
using Printf
using Plots
using FFTW

font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

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

function fps_sine(nx,ny,dx,dy,f)

    #kx = Array{Float64}(undef,nx)
    #ky = Array{Float64}(undef,ny)

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
        ue[i,j] = sin(2.0*pi*x[i])*sin(2.0*pi*y[j]) +
               c1*sin(km*2.0*pi*x[i])*sin(km*2.0*pi*y[j])

        f[i,j] = c2*sin(2.0*pi*x[i])*sin(2.0*pi*y[j]) +
                 c2*sin(km*2.0*pi*x[i])*sin(km*2.0*pi*y[j])

        un[i,j] = 0.0
    end
end

val, t, bytes, gctime, memallocs = @timed begin

un[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,f)

end

uerror = zeros(nx+1, ny+1)
rms_error = 0.0

uerror = un - ue

rms_error = compute_l2norm(nx, ny, uerror)
max_error = maximum(abs.(uerror))

println("Error details:");
println("L-2 Norm = ", rms_error);
println("Maximum Norm = ", max_error);
print("CPU Time = ", t);

p1 = contour(x, y, transpose(ue), fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Exact")
p2 = contour(x, y, transpose(un), fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Fast Sine transform")
p3 = plot(p1,p2, size = (1000, 400))
savefig(p3,"fst_contour.pdf")
