#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Single species value function generation. Will be used for the boundry condition in 2d.
The value functions solved by using the HJB equation for 1D problems. The value function 
for the prey and predator terms are computed separately.


"""

from __future__ import print_function
from fenics import *
import matplotlib.pyplot as plt
import numpy as np
from softfunctions import *

def hjb_prey1d(steps, grid, r, K, PrHer, sigmax, dc, plot,tol): 
    """
    Instruction for both functions:
    Input:
        steps - number of timesteps in the backstepping
        grid - number of grid points
        r - intrinsic rowth rate of the prey
        K - Carrying capacity
        PrHer - price of one unit of prey
        sigmax - sigma for the diffusivity term
        dc - discount rate
        plot - 1 if user wants to plot the value function (only used while programming)
        tol - tolerance

    Output:
        Vs - Value function for the prey

    """

    print('Computing the value function for the boundary conditions.')
    print('1D HJB-Predator term')
    print('')
    print('')


    max_biomass=3
    # definig mesh
    mesh = IntervalMesh(grid, 0, max_biomass)

    dt = -0.001# time step size . Note: This is a backward iteration; negative dt
    num_steps = steps   # number of time steps
    T = abs(num_steps * dt)

    # definig Function space on this mesh using Lagrange polynoimals of degree 2.
    V = FunctionSpace(mesh, "CG", 2)
    W = VectorFunctionSpace(mesh, "CG", 2)
    Vold = interpolate(Constant(0), V)

    # definign boundary values
    u0 = Expression("0",degree=1)

    # this functions checks whether the input x is on the boundary or not.
    def DirichletBoundary(x, on_boundary):
        return on_boundary & ((x[0]<1e-6))

    # Enforcing u = u0 at x = 0
    bc = DirichletBC(V, u0, DirichletBoundary)

    #prey constants
    r=r
    K=K
    PrHer = PrHer
    sigmax = sigmax
    dc=dc
    fx=interpolate(Constant(0), V)#0.2

   #upper bound on haravestign
    Fmax=5
    
    dVdx = project(Vold.dx(0),V)

    fx = Expression( 'pow(PrHer,2) / 4 / (x[0] * pow(dVdx,2) + tiny )', degree = 2,domain = mesh,dVdx=project(Vold.dx(0),V),PrHer=PrHer,tiny=pow(PrHer,2)/4/Fmax)
    # Define arbitrary drift dynamics 
    f = Expression( ( 'r*x[0]*(1-x[0]/K) - fx*x[0]'), degree=1,domain=mesh, r=r, K=K, fx=fx)
    #The diffusitivity term is written as a tensor matrix
    g = Expression( ('sigmax*sigmax*x[0]*x[0]'), degree=1,domain=mesh, sigmax=sigmax) 
    h = Expression( 'PrHer*sqrt(fx*x[0])',degree=1,domain=mesh, PrHer=PrHer, fx=fx)
    ad_dif = Expression (('sigmax*sigmax*x[0]'),degree = 2, domain = mesh, sigmax = sigmax)

    # Setting up the variational problem
    # Define functions for HJB
    Val = TrialFunction(V) # this is unkown value at time t 
    Vs = Function(V)       # This is for storing the solution 
    u = TestFunction(V)    # This is an arbitrary test function


    hjb_a = Val/dt*u*dx + dot(Val.dx(0),f-ad_dif)*u*dx - 0.5*dot(g*grad(Val),grad(u)) * dx - Val*dc*u*dx
    hjb_L = (Vold/dt - h)*u*dx #linear

    tol = tol
    nerr=1

   
    # solving the variational problem.
    #for n in range (num_steps):
    while nerr>tol:
        #update current time
        solve( hjb_a == hjb_L, Vs, bc)
        fx.dVdx = project(softpos(Vs.dx(0)))

        nerr =  norm(project(Vold-Vs,V))
        print(nerr)
        Vold.assign(Vs)


    print('')
    print('')
    print('Successfully generated the value function for the prey.')
    if plot==1:
        draw(vx,grid,max_biomass) #have to compute vertex values so this works, seems strange
   
    return(Vs)

def hjb_pred1d(steps, grid, dr, PrCod, sigmay, dc, plot,tol):
    """
    Instruction for both functions:
    Input:
        steps - number of timesteps in the backstepping
        grid - number of grid points
        dr - death rate of the predator
        PrCod - price of one unit of predator
        sigmay - sigma for the diffusivity term
        dc - discount rate
        plot - 1 if user wants to plot the value function (only used while programming)

    Output:
        Vs - Value function for the prey

    """
    

    print('Computing the value function for the boundary conditions.')
    print('1D HJB-Predator term')
    print('')
    print('')



    max_biomass=3
    # definig mesh
    mesh = IntervalMesh(grid, 0, max_biomass)

    dt = -0.001# time step size . Note: This is a backward iteration; negative dt
    num_steps = steps   # number of time steps
    T = abs(num_steps * dt)

    # defining Function space 
    V = FunctionSpace(mesh, "CG", 2)
    W = VectorFunctionSpace(mesh, "CG", 2)
    Vold = interpolate(Constant(0), V)

    # definign boundary values
    u0 = Expression("0",degree=1)
 
    # this functions checks whether the input x is on the boundary or not.
    def DirichletBoundary(x, on_boundary):
        return on_boundary & ((x[0]<1e-6))

    # Enforcing u = u0 at x = 0
    bc = DirichletBC(V, u0, DirichletBoundary)

    #predator constants
    dr=dr #death rate
    PrCod = PrCod
    sigmay = sigmay
    dc=dc
    fy=interpolate(Constant(0), V)


    #Upper bound on harvesting
    Fmax = 5

    dVdy = project(Vold.dx(0),V)
    
    fy = Expression( 'pow(PrCod,2) / 4 / (x[0] * pow(dVdy,2) + tiny )', degree = 2,domain = mesh,dVdy=project(Vold.dx(0),V),PrCod=PrCod,tiny=pow(PrCod,2)/4/Fmax)  
    # Define arbitrary drift dynamics 
    f = Expression( ( '0- dr * x[0] - fy * x[0]'), degree=1,domain=mesh, dr=dr, fy=fy)
    #The diffusitivity term is written as a tensor matrix
    g = Expression( ('sigmay*sigmay*x[0]*x[0]'), degree=1,domain=mesh, sigmay=sigmay) 
    h = Expression( 'PrCod * sqrt(fy*x[0])',degree=1,domain=mesh, PrCod=PrCod, fy=fy)
    ad_dif = Expression (('sigmay*sigmay*x[0]'),degree = 2, domain = mesh, sigmay = sigmay)


    # Setting up the variational problem
    # Define functions for HJB
    Val = TrialFunction(V) # this is unkown value at time t 
    Vs = Function(V)       # This is for storing the solution 
    u = TestFunction(V)    # This is an arbitrary test function


    hjb_a = Val/dt*u*dx + dot(Val.dx(0),f-ad_dif)*u*dx - 0.5*dot(g*grad(Val),grad(u)) * dx - Val*dc*u*dx 
    hjb_L = (Vold/dt - h)*u*dx #linear

    tol = tol
    nerr=1


    # solving the variational problem.

    while nerr>tol:
        #update current time
        solve( hjb_a == hjb_L, Vs, bc)
        fy.dVdx = project(softpos(Vs.dx(0)))

        nerr =  norm(project(Vold-Vs,V))
        print(nerr)
        Vold.assign(Vs)


    print('')
    print('')
    print('Successfully generated the value function for the predator.')
    if plot==1:
        draw(vx,grid,max_biomass) 

    return(Vs)


def draw(bleh,grid,max_biomass):
    """
    Simple plot function, input is the vertex value and final time.
    Created so the code looks muche better.

    Input:
        bleh - vertex values of the value function
        grid - number of grid points
        max_biomass - maxiumum biomass
    Output:
        Computes a plot of the value function
    """
    biomass=np.arange(0,max_biomass+(max_biomass/grid),max_biomass/grid)
    plt.plot(biomass,bleh)
    plt.xlabel('Biomass x')
    plt.ylabel('Value $V_0(x)$')
    plt.axis([0, max_biomass, 0, max(bleh)])
    plt.show()



#Function for generating the mesh in the 2d problem
def mesh_generator(lx,ly, Nx,Ny):
    m = UnitSquareMesh(Nx, Ny)
    xy = m.coordinates()

    xy=pow(xy,2)
    #Scale
    xy[:,0] = xy[:,0]*lx
    xy[:,1] = xy[:,1]*ly

    m.coordinates()[:] = xy

    return m


"""
##############################################################
# below are parameters that can be used for testing
##############################################################

nx=128
dc=0.05
r=1
K=1
Cmax = 0.5
beta = 3
PrHer = 20
sigmax = 1/sqrt(2)
epsilon=0.6
dr=0.2 #death rate
PrCod = 50
sigmay = 1/sqrt(2)
plotting=0

a=hjb_prey1d(10,nx,r,K,PrHer,sigmax,dc,plotting)
b=hjb_pred1d(10,nx,dr,PrCod,sigmay,dc,plotting)
f = Expression( ( 'a(x[0])'), degree=1,domain=mesh)
#vx=a.compute_vertex_values()
#print(max(vx))
print(a(V))
#Þetta virkar sem function!
#draw(vx,grid,max_biomass) #have to compute vertex values so this works, seems strange
#draw(vx,grid,max_biomass)
#import matplotlib.pyplot
# plotting solution
#matplotlib.pyplot.plot(u)
#matplotlib.pyplot.show()
#plot(u, interactive = True)

"""
