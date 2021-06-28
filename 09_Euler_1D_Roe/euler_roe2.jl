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
        x[i] = -0.5*dx + dx*(i)
    end
	plot(x)
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
			for i = 1:nx for m = 1:3
            	q[i,m,ri] = qn[i,m]
			end end
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

	f = Array{Float64}(undef,nx,3)

	#compute fluxes
	for i = 1:nx
		p = (gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1])
		f[i,1] = q[i,2]
		f[i,2] = q[i,2]*q[i,2]/q[i,1] + p
		f[i,3] = q[i,2]*q[i,3]/q[i,1] + p*q[i,2]/q[i,1]
	end

	# max wave speed
	for i = 1:nx
		a = sqrt((gamma*((gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1]))/q[i,1]))
		alp[i]=abs(q[i,2]/q[i,1]) + a
	end

	# max wave speed within local stencil
	for i = 3:nx-3
		a = sqrt((gamma*((gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1]))/q[i,1]))
		alpha[i]=max(alp[i-2],alp[i-1],alp[i],alp[i+1],alp[i+2])
	end
	alpha[1] = alp[1]
	alpha[2] = alp[2]
	alpha[nx-1] = alp[nx-1]
	alpha[nx-2] = alp[nx-2]


	#flux splitting
	for i = 1:nx for m=1:3
		fp[i,m] = 0.5*(f[i,m]+alpha[i]*q[i,m])
		fm[i,m] = 0.5*(f[i,m]-alpha[i]*q[i,m])
	end end


	cwp(nx,gamma,q,fp,fm,fpw,fmw)

	# WENO Reconstruction
	qL = wenoL(nx,q)
    qR = wenoR(nx,q)

	# Computing fluxes
	fluxes(nx,gamma,qL,fL)
	fluxes(nx,gamma,qR,fR)

	# compute Riemann solver using Gudanov scheme
	# gudanov(nx,gamma,qL,qR,f,fL,fR)

	# compute Riemann solver using HLLC scheme
	#hllc(nx,gamma,qL,qR,f,fL,fR)
	roe(nx,gamma,qL,qR,f,fL,fR)

	# RHS
	for i = 1:nx for m = 1:3
		r[i,m] = -(f[i+1,m] - f[i,m])/dx
	end end
end

#-----------------------------------------------------------------------------#
# Riemann solver: Roe's approximate Riemann solver
#-----------------------------------------------------------------------------#
function cwp(nx,gamma,u,fp,fm,fpw,fmw)
	dd = Array{Float64}(undef,3)
	dF = Array{Float64}(undef,3)
	V = Array{Float64}(undef,3)
	gm = gamma-1.0

	for i = 1:nx
		#Left and right states:
		rhLL = uL[i,1]
		uuLL = uL[i,2]/rhLL
		eeLL = uL[i,3]/rhLL
	    ppLL = gm*(eeLL*rhLL - 0.5*rhLL*(uuLL*uuLL))
	    hhLL = eeLL + ppLL/rhLL

		rhRR = uR[i,1]
		uuRR = uR[i,2]/rhRR
		eeRR = uR[i,3]/rhRR
	    ppRR = gm*(eeRR*rhRR - 0.5*rhRR*(uuRR*uuRR))
	    hhRR = eeRR + ppRR/rhRR

		alpha = 1.0/(sqrt(abs(rhLL)) + sqrt(abs(rhRR)))

		uu = (sqrt(abs(rhLL))*uuLL + sqrt(abs(rhRR))*uuRR)*alpha
		hh = (sqrt(abs(rhLL))*hhLL + sqrt(abs(rhRR))*hhRR)*alpha
		aa = sqrt(abs(gm*(hh-0.5*uu*uu)))

		D11 = abs(uu)
		D22 = abs(uu + aa)
		D33 = abs(uu - aa)

		beta = 0.5/(aa*aa)
		phi2 = 0.5*gm*uu*uu

		#Right eigenvector matrix
		R11, R21, R31 = 1.0, uu, phi2/gm
		R12, R22, R32 = beta, beta*(uu + aa), beta*(hh + uu*aa)
		R13, R23, R33 = beta, beta*(uu - aa), beta*(hh - uu*aa)

		#Left eigenvector matrix
		L11, L12, L13 = 1.0-phi2/(aa*aa), gm*uu/(aa*aa), -gm/(aa*aa)
		L21, L22, L23 = phi2 - uu*aa, aa - gm*uu, gm
		L31, L32, L33 = phi2 + uu*aa, -aa - gm*uu, gm

		for m = 1:3
			V[m] = 0.5*(uR[i,m]-uL[i,m])
		end

		dd[1] = D11*(L11*V[1] + L12*V[2] + L13*V[3])
		dd[2] = D22*(L21*V[1] + L22*V[2] + L23*V[3])
		dd[3] = D33*(L31*V[1] + L32*V[2] + L33*V[3])

		dF[1] = R11*dd[1] + R12*dd[2] + R13*dd[3]
		dF[2] = R21*dd[1] + R22*dd[2] + R23*dd[3]
		dF[3] = R31*dd[1] + R32*dd[2] + R33*dd[3]

		for m = 1:3
			f[i,m] = 0.5*(fR[i,m]+fL[i,m]) - dF[m]
		end
	end
end

#-----------------------------------------------------------------------------#
# Riemann solver: Roe's approximate Riemann solver
#-----------------------------------------------------------------------------#
function roe(nx,gamma,uL,uR,f,fL,fR)
	dd = Array{Float64}(undef,3)
	dF = Array{Float64}(undef,3)
	V = Array{Float64}(undef,3)
	gm = gamma-1.0

	for i = 1:nx+1
		#Left and right states:
		rhLL = uL[i,1]
		uuLL = uL[i,2]/rhLL
		eeLL = uL[i,3]/rhLL
	    ppLL = gm*(eeLL*rhLL - 0.5*rhLL*(uuLL*uuLL))
	    hhLL = eeLL + ppLL/rhLL

		rhRR = uR[i,1]
		uuRR = uR[i,2]/rhRR
		eeRR = uR[i,3]/rhRR
	    ppRR = gm*(eeRR*rhRR - 0.5*rhRR*(uuRR*uuRR))
	    hhRR = eeRR + ppRR/rhRR

		alpha = 1.0/(sqrt(abs(rhLL)) + sqrt(abs(rhRR)))

		uu = (sqrt(abs(rhLL))*uuLL + sqrt(abs(rhRR))*uuRR)*alpha
		hh = (sqrt(abs(rhLL))*hhLL + sqrt(abs(rhRR))*hhRR)*alpha
		aa = sqrt(abs(gm*(hh-0.5*uu*uu)))

		D11 = abs(uu)
		D22 = abs(uu + aa)
		D33 = abs(uu - aa)

		beta = 0.5/(aa*aa)
		phi2 = 0.5*gm*uu*uu

		#Right eigenvector matrix
		R11, R21, R31 = 1.0, uu, phi2/gm
		R12, R22, R32 = beta, beta*(uu + aa), beta*(hh + uu*aa)
		R13, R23, R33 = beta, beta*(uu - aa), beta*(hh - uu*aa)

		#Left eigenvector matrix
		L11, L12, L13 = 1.0-phi2/(aa*aa), gm*uu/(aa*aa), -gm/(aa*aa)
		L21, L22, L23 = phi2 - uu*aa, aa - gm*uu, gm
		L31, L32, L33 = phi2 + uu*aa, -aa - gm*uu, gm

		for m = 1:3
			V[m] = 0.5*(uR[i,m]-uL[i,m])
		end

		dd[1] = D11*(L11*V[1] + L12*V[2] + L13*V[3])
		dd[2] = D22*(L21*V[1] + L22*V[2] + L23*V[3])
		dd[3] = D33*(L31*V[1] + L32*V[2] + L33*V[3])

		dF[1] = R11*dd[1] + R12*dd[2] + R13*dd[3]
		dF[2] = R21*dd[1] + R22*dd[2] + R23*dd[3]
		dF[3] = R31*dd[1] + R32*dd[2] + R33*dd[3]

		for m = 1:3
			f[i,m] = 0.5*(fR[i,m]+fL[i,m]) - dF[m]
		end
	end
end

#-----------------------------------------------------------------------------#
# Riemann solver: Rusanov
#-----------------------------------------------------------------------------#
function gudanov(nx,gamma,qL,qR,f,fL,fR)

	ps = Array{Float64}(undef,nx+1)

	#wavespeed(nx,gamma,q,ps)
	wavespeed(nx,gamma,qL,qR,ps)
	# Interface fluxes (Rusanov)
	for i = 1:nx+1 for m = 1:3
		f[i,m] = 0.5*(fR[i,m]+fL[i,m]) - 0.5*ps[i]*(qR[i,m]-qL[i,m])
	end end
end

#-----------------------------------------------------------------------------#
# compute wave speed
#-----------------------------------------------------------------------------#
function wavespeed2(nx,gamma,q,ps)
	rad = Array{Float64}(undef,nx)
	# spectral radius of Jacobian
	for i = 1:nx
		a = sqrt((gamma*((gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1]))/q[i,1]))
		l1 = abs(q[i,2]/q[i,1])
		l2 = abs(q[i,2]/q[i,1] + a)
		l3 = abs(q[i,2]/q[i,1] - a)
		rad[i] = max(l1,l2,l3)
	end

	# propagation speed
	for i = 2:nx
		ps[i] = max(rad[i-1], rad[i])
	end
	ps[1] = ps[2]
	ps[nx+1] = ps[nx]
end

function wavespeed(nx,gamma,uL,uR,ps)
	rad = Array{Float64}(undef,nx)
	# spectral radius of Jacobian
	gm = gamma-1.0
	for i = 1:nx+1
		#Left and right states:
		rhLL = uL[i,1]
		uuLL = uL[i,2]/rhLL
		eeLL = uL[i,3]/rhLL
	    ppLL = gm*(eeLL*rhLL - 0.5*rhLL*(uuLL*uuLL))
	    hhLL = eeLL + ppLL/rhLL


		rhRR = uR[i,1]
		uuRR = uR[i,2]/rhRR
		eeRR = uR[i,3]/rhRR
	    ppRR = gm*(eeRR*rhRR - 0.5*rhRR*(uuRR*uuRR))
	    hhRR = eeRR + ppRR/rhRR

		alpha = 1.0/(sqrt(abs(rhLL)) + sqrt(abs(rhRR)))

		uu = (sqrt(abs(rhLL))*uuLL + sqrt(abs(rhRR))*uuRR)*alpha
		hh = (sqrt(abs(rhLL))*hhLL + sqrt(abs(rhRR))*hhRR)*alpha
		aa = sqrt(abs(gm*(hh-0.5*uu*uu)))

#		ps[i] = aa + abs(uu)
		ps[i] = abs(aa + uu)
	end
end
#-----------------------------------------------------------------------------#
# WENO reconstruction for upwind direction (positive; left to right)
# u(i): solution values at finite difference grid nodes i = 1,...,N
# f(j): reconstructed values at nodes j = i-1/2; j = 1,...,N+1
#-----------------------------------------------------------------------------#
function wenoL(n,u)
	f = Array{Float64}(undef,n+1,3)

	for m = 1:3
	    i = 0
	    v1 = u[i+3,m]
	    v2 = u[i+2,m]
	    v3 = u[i+1,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)

	    i = 1
	    v1 = u[i+1,m]
	    v2 = u[i,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)

	    i = 2
	    v1 = u[i-1,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)

	    for i = 3:n-2
	        v1 = u[i-2,m]
	        v2 = u[i-1,m]
	        v3 = u[i,m]
	        v4 = u[i+1,m]
	        v5 = u[i+2,m]
	        f[i+1,m] = wcL(v1,v2,v3,v4,v5)
	    end

	    i = n-1
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+1,m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)

	    i = n
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i,m]
	    v5 = u[i-1,m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)
	end
	return f
end

#-----------------------------------------------------------------------------#
# WENO reconstruction for downwind direction (negative; right to left)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
#-----------------------------------------------------------------------------#
function wenoR(n,u)
	f = Array{Float64}(undef,n+1,3)

	for m = 1:3
	    i = 1
	    v1 = u[i+1,m]
	    v2 = u[i,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    i = 2
	    v1 = u[i-1,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    for i = 3:n-2
	        v1 = u[i-2,m]
	        v2 = u[i-1,m]
	        v3 = u[i,m]
	        v4 = u[i+1,m]
	        v5 = u[i+2,m]
	        f[i,m] = wcR(v1,v2,v3,v4,v5)
	    end

	    i = n-1
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+1,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    i = n
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i,m]
	    v5 = u[i-1,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    i = n+1
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i-1,m]
	    v4 = u[i-2,m]
	    v5 = u[i-3,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)
	end
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
nx = 256
ns = 20
dt = 0.0001
tm = 0.20

dx = 1.0/nx
nt = Int64(tm/dt)
ds = tm/ns

#q = Array{Float64}(zero, nx,3,ns+1)
q = zeros(Float64,nx,3,ns+1)
numerical(nx,ns,nt,dx,dt,q)

x = Array(0.5*dx:dx:1.0-0.5*dx)

solution_d = open("solution_d.txt", "w")
solution_v = open("solution_v.txt", "w")
solution_e = open("solution_e.txt", "w")
for i = 1:nx
    write(solution_d, string(x[i]), " ",)
	write(solution_v, string(x[i]), " ",)
	write(solution_e, string(x[i]), " ",)
    for n = 1:ns+1
        write(solution_d, string(q[i,1,n]), " ")
		write(solution_v, string(q[i,2,n]), " ")
		write(solution_e, string(q[i,3,n]), " ")
    end
    write(solution_d, "\n",)
	write(solution_v, "\n",)
	write(solution_e, "\n",)
end
close(solution_d)
close(solution_v)
close(solution_e)
