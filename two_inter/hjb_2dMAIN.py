 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Computing the optimal harvesting policy for the stochastic Rosenzweig-MacArthur model using the HJB equation.

"""

from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from boundary_value import *
import time 

name='Base Model'
maxX=0.66 #optimal profit state of the system
maxY=0.01

start_time = time.time()

## Duration of simulation and time steps
dt = -1000.01 #time step size . Note: This is a backward iteration; negative dt

# Create mesh and define function space
#Mesh is extracted from the boundary_value.py file
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
tol=1e-6 #Convergnece tolerance

#Prey
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

#Generating the boundary condition
#The boundary values are generated by calling a function in the file boundary_value.py
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

#upper bound on the harvesti
Fmax = 3

fx = Expression( 'pow(PrHer,2) / 4 / (x[0] * pow(0.5*(dVdx+abs(dVdx)),2) + tiny )', degree = 2,domain = mesh,dVdx=project((Vold.dx(0)),V),PrHer=PrHer,tiny=pow(PrHer,2)/4/Fmax)
fy = Expression( 'pow(PrCod,2) / 4 / (x[1] * pow(0.5*(dVdy+abs(dVdy)),2) + tiny )', degree = 2,domain = mesh,dVdy=project((Vold.dx(1)),V),PrCod=PrCod,tiny=pow(PrCod,2)/4/Fmax)

# Functional response (interaction) 
Feed = Expression( 'Cmax*beta*x[0]*x[1]/(beta*x[0] + Cmax)' , degree = 2,domain = mesh,Cmax = Cmax, beta=beta)

# Define drift dynamics 
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
#imposing the boundary conditions
bc = DirichletBC(V, bv, boundary)

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

#solving the HJB equation until it converges
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


#%%
def myplot(f,title="",xlabel="Biomass Prey",ylabel="Biomass Predator"):
    plt.colorbar(plot(project(f,V)))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.axis([0, 2, 0, 2])
    plt.show()

#plotting the value function and the optimal harvesting policy
myplot(Vs,'Base model')
myplot(project(fx,V),'Prey')
myplot(project(fy,V),'Predator')



import pylab
T = 100
nsteps = 100000
#Noise Generator
rBM = lambda tv : pylab.cumsum( pylab.randn(tv.size) * pylab.sqrt(pylab.diff( pylab.append(0,tv ))))
tv = pylab.linspace(0,T,nsteps+1)
dt = tv[2]-tv[1]
Bx = rBM(tv)
By = rBM(tv)

#Define functions that are used in the population dynamics simulation
C = lambda x,y : Cmax*beta*x*y/(beta*x+Cmax)
fx1 = lambda x,y : r*x*(1-x/K) - C(x,y) - fx(x,y)*x 
fy1 = lambda x,y : epsilon*C(x,y) - dr*y - fy(x,y)*y
gx = lambda x : sigmax * x
gy = lambda y : sigmay * y

#profit functions for the prey and predator
profitX = lambda x,y : PrHer * sqrt(fx(x,y)*x)
profitY = lambda x,y : PrCod * sqrt(fy(x,y)*y)

#Catch rate functions for the prey and predator
CX = lambda x,y : fx(x,y)*x
CY = lambda x,y : fy(x,y)*y

#Allocating memory
X = pylab.zeros( tv.size )
Y = pylab.zeros( tv.size )
X1 = pylab.zeros( tv.size )
Y1 = pylab.zeros( tv.size )
H = pylab.zeros( tv.size )
cumH = pylab.zeros( tv.size )
arbH = pylab.zeros( tv.size )
CatchX =pylab.zeros( tv.size )
CatchY =pylab.zeros( tv.size )
total_catch= pylab.zeros( tv.size )
maxprof=pylab.zeros( tv.size )




X[0] = 0.3
Y[0] = 0.4
X1[0] = 0.3
Y1[0] = 0.4
cumH[0]=0
arbH[0]=0

H[0]=0
real=1 #500 simulations or it is set to one when computing only one 
X_bin=[0]*real
Y_bin= [0]*real
profit_bin=[0]*real

#Here the simulations take place
#The profit
for j in range(real):
    print(j)
    Bx = rBM(tv)
    By = rBM(tv)
    for i in range(nsteps):
        #The population dynamics for the prey (X) and predator (Y)
        X[i+1] =  X[i] + fx1(X[i],Y[i]) * dt + gx(X[i]) * (Bx[i+1]-Bx[i]) 
        Y[i+1] =  Y[i] + fy1(X[i],Y[i]) * dt + gy(Y[i]) * (By[i+1]-By[i])
        #instantaenous profit in the time step
        H[i] = (profitX(X[i], Y[i])+ profitY(X[i], Y[i]))*dt 
        #instantaenous profit is added to the cumulative profit
        cumH[i+1] = cumH[i]+H[i]

        #maximum sustainable profit that could be attained in each time step
        maxprof[i]=(profitX(maxX, maxY)+ profitY(maxX, maxY))*dt
        #The maximum sustianable cumulative profit that is possible to attain
        arbH[i+1] = arbH[i]+(profitX(maxX, maxY)+ profitY(maxX, maxY))*dt
        #Catch rate for the prey and predator also the total catch
        CatchX[i] = CX(X[i], Y[i])
        CatchY[i] = CY(X[i], Y[i])
        total_catch[i]=CX(X[i], Y[i])+CY(X[i], Y[i])
    
    #Extracting the final state of the prey and predator for each simulation for the histograms   
    X_bin[j]=X[-1]
    Y_bin[j]=Y[-1]
    #Extracting the final true profit from the simulation for the histograms   
    profit_bin[j]=cumH[-1]

#Computing the histograms that say if species is eradicated in the simulations
binwidth=1
a=plt.hist(X_bin, bins=np.arange(-0.99995, 2, binwidth))
b=a[0]
plt.xlabel('Biomass')
plt.ylabel('Count') 
plt.title(r'$Prey\ biomass:\ <=0 =%.3f,\ >0=%.3f$'%(b[0]/(b[0]+b[1]),(b[1]/(b[0]+b[1])))) 
plt.show()

a=plt.hist(Y_bin, bins=np.arange(-0.99995, 2, binwidth))
b=a[0]
plt.xlabel('Biomass')
plt.ylabel('Count')
plt.title(r'$Predator \ biomass:\ <=0 =%.3f,\ >0=%.3f$'%(b[0]/(b[0]+b[1]),(b[1]/(b[0]+b[1]))))
plt.show()

#Histogram that describes the true profit from every simulation
n, bins = np.histogram(profit_bin)
mids = 0.5*(bins[1:] + bins[:-1])
mean = np.average(mids, weights=n)
var = np.average((mids - mean)**2, weights=n)
plt.hist(profit_bin, bins='auto')
plt.xlabel('Profit')
plt.ylabel('Count')
plt.title(r'$True\ Profit\ histogram\,\ \mu=%.3f,\ \sigma=%.3f$'%(mean,np.sqrt(var)))
plt.show()


#Removing the last elements from the instantenous profit and the max sustaianble profit at each time step as they are equal to zero.
#A new time axis is constructed that is equal in length for these vectors
H=H[:-1]
maxprof=maxprof[:-1]
b=np.arange(0,(tv[2]-tv[1])*(nsteps),tv[2]-tv[1])

#The time series that staisfies the vectors besides the instantenous profit and the max sustaianble profit at each time step
a=np.arange(0,(tv[2]-tv[1])*(nsteps+1),tv[2]-tv[1])



##plotting the catch
ax3 = plt.subplot2grid((2,2), (0,0))
plt.xlabel('Time')
plt.ylabel(r'$f_xx$')
plt.plot(a, CatchX)

ax4 = plt.subplot2grid((2,2), (0,1))
plt.xlabel('Time')
plt.ylabel(r'$f_yy$')
plt.plot(a, CatchY)

ax5 = plt.subplot2grid((2,2), (1,0), colspan=2)
plt.xlabel('Time')
plt.ylabel(r'$f_xx+f_yy$')
plt.plot(a, total_catch)
plt.show()

#Computing the time series population dynamics
plt.plot(a,X,label='Prey')
plt.plot(a,Y,label='Predator')
plt.xlabel('Time')
plt.ylabel('Biomass')
plt.legend(loc='upper left')
plt.show()

#Generating the phase plane plot
xmax=max(X)+0.1 #this ensures that the axis for the plot show the whole trajectory
plt.plot(X,Y)
plt.axis([0, xmax, 0, 1])
plt.xlabel('Biomass Prey')  
plt.ylabel('Biomass Predator')

C = lambda x,y : Cmax*beta*x*y/(beta*x+Cmax)
fx1 = lambda x,y : r*x*(1-x/K) - C(x,y)
fy1 = lambda x,y : epsilon*C(x,y) - dr*y

grid=np.arange(0,1,0.01)
l_list=len(grid)

#plotting the sustaianble zone for the determenistic model
for i in range (l_list):
    for j in range(l_list):
        if fx1(grid[i],grid[j]) >0 and fy1(grid[i],grid[j]) >0:
            plt.plot(grid[i],grid[j],marker='+',color='tab:red')

plt.plot(maxX, maxY, color='lime', marker='*', linewidth=2) #plotting the optimal harvest state
plt.show()


# plotting the profit
ax1 = plt.subplot2grid((1,2), (0,0))
plt.xlabel('Time')
plt.ylabel('Profit')
plt.title('Harvest Profit')
plt.plot(b, H, label='True Profit')
plt.plot(b, maxprof, label='Max Profit')
plt.legend(loc='upper left')

ax2 = plt.subplot2grid((1,2), (0,1))
plt.xlabel('Time')
plt.ylabel('Profit')
plt.title('Cumulative Profit')
plt.plot(a, cumH,label='True Profit')
plt.plot(a,arbH, label='Max Profit')
plt.legend(loc='upper left')
plt.show()

#printing the max profit and the true profit
print('MaxProfit= ',arbH[nsteps])
print('TrueProfit= ',cumH[nsteps])
