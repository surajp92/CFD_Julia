using CPUTime
using Printf
using Plots
font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 5th-order Compact WENO scheme for spatial terms
#-----------------------------------------------------------------------------#
function numerical(nx,ns,nt,dx,dt,q)
    x = Array{Float64}(undef, nx)
    qn = Array{Float64}(undef, nx,3) # numerical solsution at every time step
    qt = Array{Float64}(undef, nx,3) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx,3)

    ri = 1 # record index
    freq = Int64(nt/ns)

    gamma = 1.4 # specific gas ratio

    # Sod's Riemann problem
    # Left side
    rhoL = 1.0
    uL = 0.0
    pL = 1.0
    # Right side
    rhoR = 0.125
    uR = 0.0
    pR = 0.1

	# nodal storage location (grid)
    for i = 1:nx
        x[i] = 0.5*dx + dx*(i-1)
    end

	xc = 0.5 # seperator location
	for i = 1:nx
	  	if (x[i] > xc)
	        rho = rhoR
			u = uR
	    	p = pR
	    else
			rho = rhoL
			u = uL
	    	p = pL
	    end

		e = p/(rho*(gamma-1.0)) + 0.5*u*u

	    #conservative variables
		qn[i,1] = rho
		qn[i,2] = rho*u
		qn[i,3] = rho*e
	end

    for i = 1:nx for m = 1:3
        q[i,m,ri] = qn[i,m] # store solution at t=0
    end	end

	# TVD RK3 for time integration
    for n = 1:nt # time step
		println(n)
        rhs(nx,dx,gamma,qn,r)

        for i = 1:nx for m = 1:3
            qt[i,m] = qn[i,m] + dt*r[i,m]
        end end

        rhs(nx,dx,gamma,qt,r)

        for i = 1:nx for m = 1:3
            qt[i,m] = 0.75*qn[i,m] + 0.25*qt[i,m] + 0.25*dt*r[i,m]
        end end

        rhs(nx,dx,gamma,qt,r)

        for i = 1:nx for m = 1:3
            qn[i,m] = (1.0/3.0)*qn[i,m] + (2.0/3.0)*qt[i,m] + (2.0/3.0)*dt*r[i,m]
        end end

        if (mod(n,freq) == 0)
            ri = ri + 1
            q[:,:,ri] = qn[:,:]
        end
    end
end

#-----------------------------------------------------------------------------#
# Calculate fluxes
#-----------------------------------------------------------------------------#
function fluxes(nx,gamma,q,f)
	for i = 1:nx+1
		p = (gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1])
		f[i,1] = q[i,2]
		f[i,2] = q[i,2]*q[i,2]/q[i,1] + p
		f[i,3] = q[i,2]*q[i,3]/q[i,1] + p*q[i,2]/q[i,1]
	end
end

#-----------------------------------------------------------------------------#
# Calculate right hand side terms of the Euler equations
#-----------------------------------------------------------------------------#

function rhs(nx,dx,gamma,q,r)
    qL = Array{Float64}(undef,nx+1,3)
    qR = Array{Float64}(undef,nx+1,3)

	fL = Array{Float64}(undef,nx+1,3)
    fR = Array{Float64}(undef,nx+1,3)

	f = Array{Float64}(undef,nx+1,3)

	rad = Array{Float64}(undef,nx)
	ps = Array{Float64}(undef,nx+1)

	# WENO Reconstruction
	for m = 1:3
  		qL[:,m] = wenoL(nx,q[:,m])
    	qR[:,m] = wenoR(nx,q[:,m])
	end

	# Computing fluxes
	fluxes(nx,gamma,qL,fL)
	fluxes(nx,gamma,qR,fL)

	for i = 1:nx
		#println(i)
		a = sqrt(gamma*((gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1]))/q[i,1])
		l1 = abs(q[i,2]/q[i,1])
		l2 = abs(q[i,2]/q[i,1] + a)
		l3 = abs(q[i,2]/q[i,1] - a)
		rad[i] = max(l1,l2,l3)
	end

	for i = 2:nx
		ps[i] = max(rad[i-1], rad[i])
	end
		ps[1] = rad[1]
		ps[nx+1] = rad[nx]

	# Interface fluxes (Rusanov)
	for i = 1:nx+1 for m = 1:3
		f[i,m] = 0.5*((fR[i,m]+fL[i,m]) - ps[i]*(qR[i,m]-qL[i,m]))
	end end

	# RHS
	for i = 1:nx for m = 1:3
		r[i,m] = -(f[i+1,m] - f[i,m])/dx
	end end
end

#-----------------------------------------------------------------------------#
# WENO reconstruction for upwind direction (positive; left to right)
# u(i): solution values at finite difference grid nodes i = 1,...,N
# f(j): reconstructed values at nodes j = i-1/2; j = 1,...,N+1
#-----------------------------------------------------------------------------#
function wenoL(n,u)
	f = Array{Float64}(undef,n+1)

    i = 1
    v1 = 4.0*u[i] - 3.0*u[i+1]
    v2 = 3.0*u[i] - 2.0*u[i+1]
    v3 = 2.0*u[i] - u[i+1]
    v4 = u[i]
    v5 = u[i+1]
    f[i] = wcL(v1,v2,v3,v4,v5)

    i = 2
    v1 = 3.0*u[i-1] - 2.0*u[i]
    v2 = 2.0*u[i-1] - u[i]
    v3 = u[i-1]
    v4 = u[i]
    v5 = u[i+1]
    f[i] = wcL(v1,v2,v3,v4,v5)

    i = 3
    v1 = 2.0*u[i-2] - u[i-1]
    v2 = u[i-2]
    v3 = u[i-1]
    v4 = u[i]
    v5 = u[i+1]
    f[i] = wcL(v1,v2,v3,v4,v5)

    for i = 4:n-1
        v1 = u[i-3]
        v2 = u[i-2]
        v3 = u[i-1]
        v4 = u[i]
        v5 = u[i+1]
        f[i] = wcL(v1,v2,v3,v4,v5)
    end

    i = n
    v1 = u[i-3]
    v2 = u[i-2]
    v3 = u[i-1]
    v4 = u[i]
    v5 = 2.0*u[i] - u[i-1]
    f[i] = wcL(v1,v2,v3,v4,v5)

    i = n+1
    v1 = u[i-3]
    v2 = u[i-2]
    v3 = u[i-1]
    v4 = 2.0*u[i-1] - u[i-2]
    v5 = 3.0*u[i-1] - 2.0*u[i-2]
    f[i] = wcL(v1,v2,v3,v4,v5)

	return f

end

#-----------------------------------------------------------------------------#
# CRWENO reconstruction for downwind direction (negative; right to left)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
#-----------------------------------------------------------------------------#
function wenoR(n,u)
	f = Array{Float64}(undef,n+1)

    i = 1
    v1 = 3.0*u[i] - 2.0*u[i+1]
    v2 = 2.0*u[i] - u[i+1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i] = wcR(v1,v2,v3,v4,v5)

    i = 2
    v1 = 2.0*u[i-1] - u[i]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i] = wcR(v1,v2,v3,v4,v5)

    for i = 3:n-2
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i]
        v4 = u[i+1]
        v5 = u[i+2]
        f[i] = wcR(v1,v2,v3,v4,v5)
    end

    i = n-1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = 2.0*u[i+1] - u[i]
    f[i] = wcR(v1,v2,v3,v4,v5)

    i = n
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = 2.0*u[i] - u[i-1]
    v5 = 3.0*u[i] - 2.0*u[i-1]
    f[i] = wcR(v1,v2,v3,v4,v5)

    i = n+1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = 2.0*u[i-1] - u[i-2]
    v4 = 3.0*u[i-1] - 2.0*u[i-2]
    v5 = 4.0*u[i-1] - 3.0*u[i-2]
    f[i] = wcR(v1,v2,v3,v4,v5)

	return f

end

#---------------------------------------------------------------------------#
#nonlinear weights for upwind direction
#---------------------------------------------------------------------------#
function wcL(v1,v2,v3,v4,v5)
    eps = 1.0e-6

    # smoothness indicators
    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^2 + 0.25*(v1-4.0*v2+3.0*v3)^2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^2 + 0.25*(v2-v4)^2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^2 + 0.25*(3.0*v3-4.0*v4+v5)^2

    # computing nonlinear weights w1,w2,w3
    c1 = 1.0e-1/((eps+s1)^2)
    c2 = 6.0e-1/((eps+s2)^2)
    c3 = 3.0e-1/((eps+s3)^2)

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    # candiate stencils
    q1 = v1/3.0 - 7.0/6.0*v2 + 11.0/6.0*v3
    q2 =-v2/6.0 + 5.0/6.0*v3 + v4/3.0
    q3 = v3/3.0 + 5.0/6.0*v4 - v5/6.0

    # reconstructed value at interface
    f = (w1*q1 + w2*q2 + w3*q3)

    return f

end

#---------------------------------------------------------------------------#
#nonlinear weights for downwind direction
#---------------------------------------------------------------------------#
function wcR(v1,v2,v3,v4,v5)
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^2 + 0.25*(v1-4.0*v2+3.0*v3)^2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^2 + 0.25*(v2-v4)^2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^2 + 0.25*(3.0*v3-4.0*v4+v5)^2

    c1 = 3.0e-1/(eps+s1)^2
    c2 = 6.0e-1/(eps+s2)^2
    c3 = 1.0e-1/(eps+s3)^2

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    # candiate stencils
    q1 =-v1/6.0      + 5.0/6.0*v2 + v3/3.0
    q2 = v2/3.0      + 5.0/6.0*v3 - v4/6.0
    q3 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0

    # reconstructed value at interface
    f = (w1*q1 + w2*q2 + w3*q3)
end

#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
nx = 512
ns = 20
dt = 0.0001
tm = 0.20

dx = 1.0/nx
nt = Int64(tm/dt)
ds = tm/ns

q = Array{Float64}(undef, nx,3,ns+1)
numerical(nx,ns,nt,dx,dt,q)

x = Array(0.5*dx:dx:1.0-0.5*dx)

qplot = q[:,1,:]
p1 = plot(x,qplot,lw = 1,
          xlabel="\$X\$", ylabel = "\$U\$",
          xlims=(minimum(x),maximum(x)),
          grid=(:none), legend=:none)

plot(p1, size = (1000, 600))
savefig("shocktube.pdf")
