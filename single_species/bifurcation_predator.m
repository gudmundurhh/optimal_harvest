%Bifurcation with varying mortality rate for the predator
clear all; clc; close all;

%Set initial conditions

P0=1;%number of predators
N0=1;%number of preys
tspan=linspace(0,800,5000);
Cmax=1;
K=1;
beta=3;
r=1;
epsilon=0.6;
dr=linspace(0,1,300);

%allocate memory
N=zeros(1,length(dr));
P=zeros(1,length(dr));
max_N=zeros(1,length(dr));
min_N=zeros(1,length(dr));
max_P=zeros(1,length(dr));
min_P=zeros(1,length(dr));

%Simulate for every mortality rate
for i=1:length(dr)
drnow=dr(i);

odefun=@(tspan,x)[r*x(1)*(1-x(1)/K)-(Cmax*beta*x(1)*x(2)/(beta*x(1)+Cmax));epsilon*(Cmax*beta*x(1)*x(2)/(beta*x(1)+Cmax))-drnow*x(2)];
[~,NP] = ode45(odefun,tspan,[N0,P0]);

%Finding the maximum and minimum rates from the last 400 points in the
%dynamics for prey and predator to see how the limit cycle behaves or if it
%is converging to equilibrium
max_N(i)=max(NP(end-400:end,1));
min_N(i)=min(NP(end-400:end,1));
max_P(i)=max(NP(end-400:end,2));
min_P(i)=min(NP(end-400:end,2));

end

%plotting the bifucation diagram
plot(dr,max_N,'b','LineWidth',2)

hold on 
plot(dr,min_P,'r','LineWidth',2)
plot(dr,min_N,'b','LineWidth',2)
axis([0,1,0,1.2])

plot(dr,max_P,'r','LineWidth',2)
xlabel('Mortality rate \mu')
ylabel('Biomass')
legend('Prey','Predator') 
title('Bifurcation Prey and Predator')

