using CPUTime
using Printf
using Plots
font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

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
function tdms(a,b,c,r,x,s,e)
    gam = Array{Float64}(undef, e)
    bet = b[s]
    x[s] = r[s]/bet

    for i = s+1:e
        gam[i] = c[i-1]/bet
        bet = b[i] - a[i]*gam[i]
        x[i] = (r[i] - a[i]*x[i-1])/bet
    end

    for i = e-1:-1:s
        x[i] = x[i] - gam[i+1]*x[i+1]
    end
end

#-----------------------------------------------------------------------------#
# Solution to tridigonal system using cyclic Thomas algorithm
#-----------------------------------------------------------------------------#
function ctdms(a,b,c,alpha,beta,r,x,s,e)
    bb = Array{Float64}(undef, e)
    u = Array{Float64}(undef, e)
    z = Array{Float64}(undef, e)
    gamma = -b[s]
    bb[s] = b[s] -gamma
    bb[e] = b[e] - alpha*beta/gamma

    for i = s+1:e-1
        bb[i] = b[i]
    end

    tdms(a,bb,c,r,x,s,e)

    u[s] = gamma
    u[e] = alpha
    for i = s+1:e-1
        u[i] = 0.0
    end

    tdms(a,bb,c,u,z,s,e)

    fact = (x[s] + beta*x[e]/gamma)/(1.0 + z[s] + beta*z[e]/gamma)

    for i = s:e
        x[i] = x[i] - fact*z[i]
    end
end

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 5th-order Compact WENO scheme for spatial terms
#-----------------------------------------------------------------------------#
function numerical(nx,ns,nt,dx,dt,u)
    x = Array{Float64}(undef, nx+1)
    un = Array{Float64}(undef, nx+1) # numerical solsution at every time step
    ut = Array{Float64}(undef, nx+1) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx)

    k = 1 # record index
    freq = Int64(nt/ns)

    for i = 1:nx+1
        x[i] = dx*(i-1)
        un[i] = sin(2.0*pi*x[i])
        u[i,k] = un[i] # store solution at t=0
    end

    for j = 2:nt+1
        rhs(nx,dx,un,r)

        for i = 1:nx
            ut[i] = un[i] + dt*r[i]
        end
        ut[nx+1] = ut[1] # periodic

        rhs(nx,dx,ut,r)

        for i = 1:nx
            ut[i] = 0.75*un[i] + 0.25*ut[i] + 0.25*dt*r[i]
        end
        ut[nx+1] = ut[1] # periodic

        rhs(nx,dx,ut,r)

        for i = 1:nx
            un[i] = (1.0/3.0)*un[i] + (2.0/3.0)*ut[i] + (2.0/3.0)*dt*r[i]
        end
        un[nx+1] = un[1] # periodic

        if (mod(j,freq) == 0)
            u[:,k] = un[:]
            k = k+1
        end
    end
end

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -u∂u/∂x
#-----------------------------------------------------------------------------#
function rhs(nx,dx,u,r)
    uL = Array{Float64}(undef, nx)
    uR = Array{Float64}(undef, nx+1)

    crwenoL(nx,u,uL)

    crwenoR(nx,u,uR)

    for i = 2:nx
        if (u[i] >= 0.0)
            r[i] = -u[i]*(uL[i] - uL[i-1])/dx
        else
            r[i] = -u[i]*(uR[i+1] - uR[i])/dx
        end
    end
    #for i = 1; periodic
    i = 1
    if (u[i] >= 0.0)
        r[i] = -u[i]*(uL[i] - uL[nx])/dx
    else
        r[i] = -u[i]*(uR[i+1] - uR[nx+1])/dx
    end
end

#-----------------------------------------------------------------------------#
# CRWENO reconstruction for upwind direction (positive; left to right)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j <== i+1/2; only use j = 1,2,...,N
#-----------------------------------------------------------------------------#
function crwenoL(n,u,f)
    a = Array{Float64}(undef, n)
    b = Array{Float64}(undef, n)
    c = Array{Float64}(undef, n)
    r = Array{Float64}(undef, n)

    i = 1
    v1 = u[n-1]
    v2 = u[n]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]

    a1,a2,a3,b1,b2,b3 = crwcL(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[n] + b2*u[i] + b3*u[i+1]

    i = 2
    v1 = u[n]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]

    a1,a2,a3,b1,b2,b3 = crwcL(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]

    for i = 3:n-1
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i]
        v4 = u[i+1]
        v5 = u[i+2]

        a1,a2,a3,b1,b2,b3 = crwcL(v1,v2,v3,v4,v5)
        a[i] = a1
        b[i] = a2
        c[i] = a3
        r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]
    end

    i = n
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[2]

    a1,a2,a3,b1,b2,b3 = crwcL(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]

    alpha = c[n]
    beta = a[1]

    ctdms(a,b,c,alpha,beta,r,f,1,n)

end

#-----------------------------------------------------------------------------#
# CRWENO reconstruction for downwind direction (negative;right to left)
# u(i): solution values at finite difference grid nodes i =1,...,N+1
# f(j): reconstructed values at nodes j <== i-1/2; only use j = 2,...,N+1
#-----------------------------------------------------------------------------#
function crwenoR(n,u,f)
    a = Array{Float64}(undef, n+1)
    b = Array{Float64}(undef, n+1)
    c = Array{Float64}(undef, n+1)
    r = Array{Float64}(undef, n+1)

    i = 2
    v1 = u[n]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]

    a1,a2,a3,b1,b2,b3 = crwcR(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]


    for i = 3:n-1
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i]
        v4 = u[i+1]
        v5 = u[i+2]

        a1,a2,a3,b1,b2,b3 = crwcR(v1,v2,v3,v4,v5)
        a[i] = a1
        b[i] = a2
        c[i] = a3
        r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]
    end

    i = n
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[2]

    a1,a2,a3,b1,b2,b3 = crwcR(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]

    i = n+1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[2]
    v5 = u[3]

    a1,a2,a3,b1,b2,b3 = crwcR(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[i-1] + b2*u[i] + b3*u[2]

    alpha = c[n+1]
    beta = a[2]

    ctdms(a,b,c,alpha,beta,r,f,2,n+1)

end

#---------------------------------------------------------------------------#
#nonlinear weights for upwind direction
#---------------------------------------------------------------------------#
function crwcL(v1,v2,v3,v4,v5)
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^2 + 0.25*(v1-4.0*v2+3.0*v3)^2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^2 + 0.25*(v2-v4)^2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^2 + 0.25*(3.0*v3-4.0*v4+v5)^2

    c1 = 2.0e-1/((eps+s1)^2)
    c2 = 5.0e-1/((eps+s2)^2)
    c3 = 3.0e-1/((eps+s3)^2)

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    a1 = (2.0*w1 + w2)/3.0
    a2 = (w1 + 2.0*w2 + 2.0*w3)/3.0
    a3 = w3/3.0

    b1 = w1/6.0
    b2 = (5.0*w1 + 5.0*w2 + w3)/6.0
    b3 = (w2 + 5.0*w3)/6.0

    return a1,a2,a3,b1,b2,b3

end

#---------------------------------------------------------------------------#
#nonlinear weights for downwind direction
#---------------------------------------------------------------------------#
function crwcR(v1,v2,v3,v4,v5)
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^2 + 0.25*(v1-4.0*v2+3.0*v3)^2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^2 + 0.25*(v2-v4)^2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^2 + 0.25*(3.0*v3-4.0*v4+v5)^2

    c1 = 3.0e-1/(eps+s1)^2
    c2 = 5.0e-1/(eps+s2)^2
    c3 = 2.0e-1/(eps+s3)^2

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    a1 = w1/3.0
    a2 = (w3 + 2.0*w2 + 2.0*w1)/3.0
    a3 = (2.0*w3 + w2)/3.0

    b1 = (w2 + 5.0*w1)/6.0
    b2 = (5.0*w3 + 5.0*w2 + w1)/6.0
    b3 = w3/6.0

    return a1,a2,a3,b1,b2,b3
end

#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
nx = 200
ns = 10
dt = 0.0001
tm = 0.25

dx = 1.0/nx
nt = Int64(tm/dt)
ds = tm/ns

u = Array{Float64}(undef, nx+1, ns+1)
numerical(nx,ns,nt,dx,dt,u)

x = Array(0:dx:1.0)

solution = open("solution_p.txt", "w")

for i = 1:nx+1
    write(solution, string(x[i]), " ",)
    for j = 1:ns
        write(solution, string(u[i,j]), " ")
    end
    write(solution, "\n",)
end

close(solution)
