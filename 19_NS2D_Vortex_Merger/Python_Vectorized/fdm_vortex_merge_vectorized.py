#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 11:08:41 2019

@author: Suraj Pawar
Two-dimensional navier-stokes solver  
Vorticity-stream function formulation
Arakawa scheme (or compact scheme or explicit) for nonlinear term
3rd order Runge-Kutta for temporal discritization
Periodic boundary conditions only

"""
import numpy as np
from numpy.random import seed
seed(1)
import pyfftw
from scipy import integrate
from scipy import linalg
import matplotlib.pyplot as plt 
import time as tm
import matplotlib.ticker as ticker

font = {'family' : 'Times New Roman',
        'size'   : 14}    
plt.rc('font', **font)
#%%
# fast poisson solver using second-order central difference scheme
def fps(nx, ny, dx, dy, f):
    epsilon = 1.0e-6
    aa = -2.0/(dx*dx) - 2.0/(dy*dy)
    bb = 2.0/(dx*dx)
    cc = 2.0/(dy*dy)
    hx = 2.0*np.pi/np.float64(nx)
    hy = 2.0*np.pi/np.float64(ny)
    
    kx = np.empty(nx)
    ky = np.empty(ny)
    
    kx[:] = hx*np.float64(np.arange(0, nx))

    ky[:] = hy*np.float64(np.arange(0, ny))
    
    kx[0] = epsilon
    ky[0] = epsilon

    kx, ky = np.meshgrid(np.cos(kx), np.cos(ky), indexing='ij')
    
    data = np.empty((nx,ny), dtype='complex128')
    data1 = np.empty((nx,ny), dtype='complex128')
    
    data[:,:] = np.vectorize(complex)(f[1:nx+1,1:ny+1],0.0)

    a = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')
    b = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')
    
    fft_object = pyfftw.FFTW(a, b, axes = (0,1), direction = 'FFTW_FORWARD')
    fft_object_inv = pyfftw.FFTW(a, b,axes = (0,1), direction = 'FFTW_BACKWARD')
    
    e = fft_object(data)
    #e = pyfftw.interfaces.scipy_fftpack.fft2(data)
    
    e[0,0] = 0.0
    
    data1[:,:] = e[:,:]/(aa + bb*kx[:,:] + cc*ky[:,:])

    ut = np.real(fft_object_inv(data1))
    
    #periodicity
    u = np.empty((nx+3,ny+3)) 
    u[1:nx+1,1:ny+1] = ut
    u[:,ny+1] = u[:,1]
    u[nx+1,:] = u[1,:]
    u[nx+1,ny+1] = u[1,1]
    return u

#%%
# set periodic boundary condition for ghost nodes. Index 0 and (n+2) are the ghost boundary locations
def bc(nx,ny,u):
    u[:,0] = u[:,ny]
    u[:,ny+2] = u[:,2]
    
    u[0,:] = u[nx,:]
    u[nx+2,:] = u[2,:]
    
    return u  
    
#%% 
# compute rhs using arakawa scheme
# computed at all physical domain points (1:nx+1,1:ny+1; all boundary points included)
# no ghost points
def rhs(nx,ny,dx,dy,re,w,s):
    aa = 1.0/(dx*dx)
    bb = 1.0/(dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0
    
    f = np.zeros((nx+3,ny+3))
    
    #Arakawa
    j1 = gg*( (w[2:nx+3,1:ny+2]-w[0:nx+1,1:ny+2])*(s[1:nx+2,2:ny+3]-s[1:nx+2,0:ny+1]) \
             -(w[1:nx+2,2:ny+3]-w[1:nx+2,0:ny+1])*(s[2:nx+3,1:ny+2]-s[0:nx+1,1:ny+2]))

    j2 = gg*( w[2:nx+3,1:ny+2]*(s[2:nx+3,2:ny+3]-s[2:nx+3,0:ny+1]) \
            - w[0:nx+1,1:ny+2]*(s[0:nx+1,2:ny+3]-s[0:nx+1,0:ny+1]) \
            - w[1:nx+2,2:ny+3]*(s[2:nx+3,2:ny+3]-s[0:nx+1,2:ny+3]) \
            + w[1:nx+2,0:ny+1]*(s[2:nx+3,0:ny+1]-s[0:nx+1,0:ny+1]))

    j3 = gg*( w[2:nx+3,2:ny+3]*(s[1:nx+2,2:ny+3]-s[2:nx+3,1:ny+2]) \
            - w[0:nx+1,0:ny+1]*(s[0:nx+1,1:ny+2]-s[1:nx+2,0:ny+1]) \
            - w[0:nx+1,2:ny+3]*(s[1:nx+2,2:ny+3]-s[0:nx+1,1:ny+2]) \
            + w[2:nx+3,0:ny+1]*(s[2:nx+3,1:ny+2]-s[1:nx+2,0:ny+1]) )

    jac = (j1+j2+j3)*hh
    
    lap = aa*(w[2:nx+3,1:ny+2]-2.0*w[1:nx+2,1:ny+2]+w[0:nx+1,1:ny+2]) \
        + bb*(w[1:nx+2,2:ny+3]-2.0*w[1:nx+2,1:ny+2]+w[1:nx+2,0:ny+1])
        
    f[1:nx+2,1:ny+2] = -jac + lap/re 
                        
    return f

#%%
# set initial condition for vortex merger problem
def vm_ic(nx,ny,x,y):
    w = np.empty((nx+3,ny+3))
    sigma = np.pi
    xc1 = np.pi-np.pi/4.0
    yc1 = np.pi
    xc2 = np.pi+np.pi/4.0
    yc2 = np.pi
    
    w[1:nx+2, 1:ny+2] = np.exp(-sigma*((x[0:nx+1, 0:ny+1]-xc1)**2 + (y[0:nx+1, 0:ny+1]-yc1)**2)) \
                        + np.exp(-sigma*((x[0:nx+1, 0:ny+1]-xc2)**2 + (y[0:nx+1, 0:ny+1]-yc2)**2))
    
    w = bc(nx,ny,w)

    return w
  
#%% 
# read input file
l1 = []
with open('input.txt') as f:
    for l in f:
        l1.append((l.strip()).split("\t"))

nd = np.int64(l1[0][0])
nt = np.int64(l1[1][0])
re = np.float64(l1[2][0])
dt = np.float64(l1[3][0])
ns = np.int64(l1[4][0])
isolver = np.int64(l1[5][0])
isc = np.int64(l1[6][0])
ich = np.int64(l1[7][0])
ipr = np.int64(l1[8][0])
ndc = np.int64(l1[9][0])

freq = int(nt/ns)

if (ich != 19):
    print("Check input.txt file")

#%% 
# assign parameters
nx = nd
ny = nd

nxc = ndc
nyc = ndc

#%%
pi = np.pi
lx = 2.0*pi
ly = 2.0*pi

dx = lx/np.float64(nx)
dy = ly/np.float64(ny)

dxc = lx/np.float64(nxc)
dyc = ly/np.float64(nyc)

ifile = 0
time = 0.0

x = np.linspace(0.0,2.0*np.pi,nx+1)
y = np.linspace(0.0,2.0*np.pi,ny+1)

x, y = np.meshgrid(x, y, indexing='ij')

#%% 
# allocate the vorticity and streamfunction arrays
w = np.empty((nx+3,ny+3)) 
s = np.empty((nx+3,ny+3))

t = np.empty((nx+3,ny+3))

r = np.empty((nx+3,ny+3))

#%%
# set the initial condition based on the problem selected
w0 = vm_ic(nx,ny,x,y)
    
w = np.copy(w0)
s = fps(nx, ny, dx, dy, -w)
s = bc(nx,ny,s)

#%%
# time integration using third-order Runge Kutta method
aa = 1.0/3.0
bb = 2.0/3.0
clock_time_init = tm.time()
for k in range(1,nt+1):
    time = time + dt
    r = rhs(nx,ny,dx,dy,re,w,s)
    
    #stage-1
    t[1:nx+2,1:ny+2] = w[1:nx+2,1:ny+2] + dt*r[1:nx+2,1:ny+2]
    
    t = bc(nx,ny,t)
    
    s = fps(nx, ny, dx, dy, -t)
    s = bc(nx,ny,s)
    
    r = rhs(nx,ny,dx,dy,re,t,s)
    
    #stage-2
    t[1:nx+2,1:ny+2] = 0.75*w[1:nx+2,1:ny+2] + 0.25*t[1:nx+2,1:ny+2] + 0.25*dt*r[1:nx+2,1:ny+2]
    
    t = bc(nx,ny,t)
    
    s = fps(nx, ny, dx, dy, -t)
    s = bc(nx,ny,s)
    
    r = rhs(nx,ny,dx,dy,re,t,s)
    
    #stage-3
    w[1:nx+2,1:ny+2] = aa*w[1:nx+2,1:ny+2] + bb*t[1:nx+2,1:ny+2] + bb*dt*r[1:nx+2,1:ny+2]
    
    w = bc(nx,ny,w)
    
    s = fps(nx, ny, dx, dy, -w)
    s = bc(nx,ny,s)
    
    if (k%freq == 0):
        #u,v = compute_velocity(nx,ny,dx,dy,s)
        #compute_stress(nx,ny,nxc,nyc,dxc,dyc,u,v,k,freq)
        #write_data(nx,ny,dx,dy,nxc,nyc,dxc,dyc,w,s,k,freq)
        print(k, " ", time)

total_clock_time = tm.time() - clock_time_init
print('Total clock time=', total_clock_time)


#%%
# contour plot for initial and final vorticity
fig, axs = plt.subplots(1,2,sharey=True,figsize=(9,5))

cs = axs[0].contourf(w0[1:nx+2,1:ny+2].T, 120, cmap = 'jet', interpolation='bilinear')
axs[0].text(0.4, -0.1, '$t = 0.0$', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top')
cs = axs[1].contourf(w[1:nx+2,1:ny+2].T, 120, cmap = 'jet', interpolation='bilinear')
axs[1].text(0.4, -0.1, '$t = '+str(dt*nt)+'$', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top')

fig.tight_layout() 

fig.subplots_adjust(bottom=0.15)

cbar_ax = fig.add_axes([0.22, -0.05, 0.6, 0.04])
fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
plt.show()

fig.savefig("field_fdm.png", bbox_inches = 'tight')
   
