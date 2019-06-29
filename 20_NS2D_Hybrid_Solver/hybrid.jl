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
    rms = 0.0
    # println(residual)
    for j = 1:ny+1 for i = 1:nx+1
        rms = rms + r[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx+1)*(ny+1)))
    return rms
end

function wavespace(nx,ny,dx,dy)
    eps = 1.0e-6

    kx = Array{Float64}(undef,nx)
    ky = Array{Float64}(undef,ny)

    k2 = Array{Float64}(undef,nx,ny)

    #wave number indexing
    hx = 2.0*pi/(nx*dx)

    for i = 1:Int64(nx/2)
        kx[i] = hx*(i-1.0)
        kx[i+Int64(nx/2)] = hx*(i-Int64(nx/2)-1)
    end
    kx[1] = eps
    ky = kx

    for i = 1:nx for j = 1:ny
        k2[i,j] = kx[i]^2 + ky[j]^2
    end end

    return k2
end
#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
function numerical(nx,ny,nt,dx,dy,dt,re,wn,ns)

    w1f = Array{Complex{Float64}}(undef,nx,ny) # temporary array during RK3 integration
    w2f = Array{Complex{Float64}}(undef,nx,ny)
    wnf = Array{Complex{Float64}}(undef,nx,ny)

    j1f = Array{Complex{Float64}}(undef,nx,ny)
    j2f = Array{Complex{Float64}}(undef,nx,ny)
    jnf = Array{Complex{Float64}}(undef,nx,ny)

    ut = Array{Float64}(undef, nx+1, ny+1)
    data = Array{Complex{Float64}}(undef,nx,ny)
    d1 = Array{Float64}(undef, nx, ny)
    d2 = Array{Float64}(undef, nx, ny)
    d3 = Array{Float64}(undef, nx, ny)

    k2 = Array{Float64}(undef, nx, ny)

    m = 1 # record index
    freq = Int64(nt/ns)

    for i = 1:nx for j = 1:ny
        data[i,j] = complex(wn[i+1,j+1],0.0)
    end end

    wnf = fft(data)
    k2 = wavespace(nx,ny,dx,dy)
    wnf[1,1] = 0.0

    alpha1, alpha2, alpha3 = 8.0/15.0, 2.0/15.0, 1.0/3.0
    gamma1, gamma2, gamma3 = 8.0/15.0, 5.0/12.0, 3.0/4.0
    rho2, rho3 = -17.0/60.0, -5.0/12.0

    for i = 1:nx for j = 1:ny
        z = 0.5*dt*k2[i,j]/re
        d1[i,j] = alpha1*z
        d2[i,j] = alpha2*z
        d3[i,j] = alpha3*z
    end end

    for k = 1:nt
        jnf = jacobian(nx,ny,dx,dy,wnf,k2)

        for i = 1:nx for j = 1:ny
            w1f[i,j] = ((1.0 - d1[i,j])/(1.0 + d1[i,j]))*wnf[i,j] +
                        (gamma1*dt*jnf[i,j])/(1.0 + d1[i,j])
        end end

        w1f[1,1] = 0.0
        j1f = jacobian(nx,ny,dx,dy,w1f,k2)

        for i = 1:nx for j = 1:ny
            w2f[i,j] = ((1.0 - d2[i,j])/(1.0 + d2[i,j]))*w1f[i,j] +
                        (rho2*dt*jnf[i,j] + gamma2*dt*j1f[i,j])/(1.0 + d2[i,j])
        end end

        w2f[1,1] = 0.0
        j2f = jacobian(nx,ny,dx,dy,w2f,k2)

        for i = 1:nx for j = 1:ny
            wnf[i,j] = ((1.0 - d3[i,j])/(1.0 + d3[i,j]))*w2f[i,j] +
                        (rho3*dt*j1f[i,j] + gamma3*dt*j2f[i,j])/(1.0 + d3[i,j])
        end end

        if (mod(k,freq) == 0)
            println(k)
            ut[1:nx,1:ny] = real(ifft(wnf))
            # periodic BC
            ut[nx+1,:] = ut[1,:]
            ut[:,ny+1] = ut[:,1]
            field_final = open(string("vm",string(m),".txt"), "w");
            for j = 1:ny+1 for i = 1:nx+1
                write(field_final, string(x[i]), " ",string(y[j]), " ", string(ut[i,j]), " \n")
            end end
            m = m+1
            close(field_final)
        end
    end

    return ut
end

#-----------------------------------------------------------------------------#
# Calculate Jacobian in fourier space
# jf = -J(w,ψ)
#-----------------------------------------------------------------------------#
function jacobian(nx,ny,dx,dy,wf,k2)

    s = Array{Float64}(undef, nx+2, ny+2)
    sf = Array{Complex{Float64}}(undef,nx,ny)
    w = Array{Float64}(undef, nx+2, ny+2)
    data = Array{Complex{Float64}}(undef,nx,ny)
    jf = Array{Complex{Float64}}(undef,nx,ny)

    # Arakawa numerical scheme for Jacobian
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0

    w[2:nx+1, 2:ny+1] = real(ifft(wf))

    # periodic BC
    w[nx+2,:] = w[2,:]
    w[:,ny+2] = w[:,2]

    # ghost points
    w[1,:] = w[nx+1,:]
    w[:,1] = w[:,ny+1]

    for i = 1:nx for j = 1:ny
        sf[i,j] = wf[i,j]/k2[i,j]
    end end

    s[2:nx+1, 2:ny+1] = real(ifft(sf))

    # periodic BC
    s[nx+2,:] = s[2,:]
    s[:,ny+2] = s[:,2]

    # ghost points
    s[1,:] = s[nx+1,:]
    s[:,1] = s[:,ny+1]

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

        data[i-1,j-1] = complex(-(j1+j2+j3)*hh, 0.0)
    end end

    jf = fft(data)

    return jf
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

un = numerical(nx,ny,nt,dx,dy,dt,re,wn,ns)

time = tf

field_final = open("field_final.txt", "w");
for j = 1:ny+1 for i = 1:nx+1
    write(field_final, string(x[i]), " ",string(y[j]), " ", string(un[i,j]), " \n")
end end

close(field_final)

p1 = contour(x, y, transpose(un), fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Numerical")
savefig(p1,"vm.pdf")
