 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A function that computes the base case again. Only used when it is compared to other cases as the original base model is not a callable function
The original can be found in the script hjb_2d.py

"""
from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from boundary_value import *
import time 
def base():


    start_time = time.time()

    ## Duration of simulation and time steps
    dt = -1000.01 #time step size 

    # Create mesh and define function space
    nx = ny = 32
    lx = ly = 2
    mesh=mesh_generator(lx,ly, nx,ny)

    V = FunctionSpace(mesh, 'CG', 1)

    # Define functions for HJB
    Val = TrialFunction(V) # this is unkown value at time t 
    Vs = Function(V)       # This is for storing the solution 
    Vold = Function(V)
    u = TestFunction(V)    # This is an arbitrary test function

    dc=0.05 #Discounting rate
    tol=1e-6
    #prey
    r=1
    K=1
    Cmax = 1
    beta = 3
    PrHer = 1e0
    sigmax = 1/sqrt(5)
    sigmax0 = 0e-6

    #predator
    epsilon=0.6
    dr=0.15 #death rate
    PrCod = 1e0
    sigmay = 1/sqrt(5)
    sigmay0 = 0e-6

    Plot_or_not = 0

    prey1d=hjb_prey1d(10, nx, r, K, PrHer, sigmax, dc, Plot_or_not,1e-3)
    pred1d=hjb_pred1d(10, nx, dr, PrCod, sigmay, dc, Plot_or_not,1e-3)

    dte = Expression('dt',degree=2,domain=mesh,dt=dt)

    class BoundaryValues(UserExpression): 
        def eval(self, values, x): 
            values[0] =  prey1d(x[0])+pred1d(x[1]) 
         
    bv = BoundaryValues()   

    # Set the terminal value function to zero
    Vs=interpolate(Constant(0), V)
    Vold=interpolate(Constant(0), V)

    dVdx = project(Vold.dx(0),V)
    dVdy = project(Vold.dx(1),V)

    Fmax = 3

    fx = Expression( 'pow(PrHer,2) / 4 / (x[0] * pow(0.5*(dVdx+abs(dVdx)),2) + tiny )', degree = 2,domain = mesh,dVdx=project((Vold.dx(0)),V),PrHer=PrHer,tiny=pow(PrHer,2)/4/Fmax)
    fy = Expression( 'pow(PrCod,2) / 4 / (x[1] * pow(0.5*(dVdy+abs(dVdy)),2) + tiny )', degree = 2,domain = mesh,dVdy=project((Vold.dx(1)),V),PrCod=PrCod,tiny=pow(PrCod,2)/4/Fmax)

    # Functional response
    Feed = Expression( 'Cmax*beta*x[0]*x[1]/(beta*x[0] + Cmax)' , degree = 2,domain = mesh,Cmax = Cmax, beta=beta)

    # Define arbitrary drift dynamics 
    f = Expression( ( 'r*x[0]*(1-x[0]/K) - fx*x[0]-Feed', 'epsilon * Feed - dr * x[1] - fy * x[1]'), degree=2,domain=mesh, r=r, K=K, Feed=Feed, epsilon=epsilon, dr=dr, fx=fx, fy=fy)

    #The diffusitivity term is written as a tensor matrix
    g = Expression( (('sigmax0 + sigmax*sigmax*x[0]*x[0]','0'), ('0','sigmay0 + sigmay*sigmay*x[1]*x[1]')), degree=2,domain=mesh, sigmax=sigmax, sigmay=sigmay,sigmax0=sigmax0,sigmay0=sigmay0)

    # Define arbitrary pay-off
    h = Expression( 'PrHer*sqrt(fx*x[0]) + PrCod * sqrt(fy*x[1])',degree=2,domain=mesh, PrHer=PrHer, fx=fx, PrCod=PrCod, fy=fy) 

    #Advection diffusion fix
    ad_dif = Expression (('sigmax*sigmax*x[0]','sigmay*sigmay*x[1]'),degree = 2, domain = mesh, sigmax = sigmax, sigmay = sigmay)

    hjb_a = Val/dte*u*dx + dot(grad(Val),f-ad_dif)*u*dx - 0.5*dot(g*grad(Val),grad(u)) * dx - Val*dc*u*dx
    hjb_L = (Vold/dte - h)*u*dx #linear


    def boundary(x, on_boundary):
        return on_boundary & ((x[0]<1e-6) | (x[1] < 1e-6 )) 

    bc = DirichletBC(V, bv, boundary)

    def myplot(f,title="",xlabel="Biomass Prey",ylabel="Biomass Predator"):
        plt.colorbar(plot(project(f,V)))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.axis([0, 2, 0, 2])
        plt.show()
        

    #%%

    Vmesh = Vold.compute_vertex_values()
    nerr=1
    ## Smoother 
    kSmooth = 0.01
    Vrough = Function(V)
    Vsmooth = TrialFunction(V)
    Vsmoothed = Function(V)
    aSmooth = (u*Vsmooth + kSmooth * inner(grad(u),grad(Vsmooth))) * dx
    bSmooth = u*Vrough * dx

    i=0
    while nerr>tol:
        i+=1
        print('Run=',i)
        Vold.assign(Vs)

        fx.dVdx = project(Vs.dx(0),V)
        fy.dVdy = project(Vs.dx(1),V)

        solve(hjb_a == hjb_L,Vs,bc)

        Vrough.assign(Vs)
        solve(aSmooth==bSmooth,Vs,bc)
        nerr =  norm(project(Vold-Vs,V))
        print(nerr)
        
        Voldmesh = Vmesh
        Vmesh = Vs.compute_vertex_values()



    elapsed_time = time.time() - start_time

    print(elapsed_time)

    return(Vs,fx,fy)





