%Optimal harvesting policy for a stochastic logistic growth equation with
%different discounting rates

clc; clear all; close all;

%Parameters
K=1; %Carrying Capacity
r=1; %Growth rate
sigma=1/sqrt(2); %Noise in population dynamics
x0 = 1; %Initial value
dr=[ 0, 0.02, 0.05, 0.10, 0.15];%Discount rates

dt=0.0001; % time steps for euler method
T=100; %time span
tvec=0:dt:T;

%Upperbound on fisheris mortality
Umax=10;

%Discretization of state space
XMax=4; %Biomass
dx=0.05; %looping over 0.05 time steps
xi=0:dx:XMax;
xc=xi(1:length(xi)-1)+0.025;

% Functions entering into the model
%Uncontrolled system advection-diffusion equation
f = @(x) r.*x.*(1-x./K);

D = @(x) 1/2*sigma.^2*x.^2;
dD = @(x) sigma.^2*x;

v = @(x) f(x)-dD(x);

G = fvade(v,D,xi,'r');


% Effect of the fishing: The "generator" d/dx
ddxl = -fvade(@(x)-1,@(x)0,xi,'r');
ddxl(1,:)=ddxl(2,:);

% simulation control
%arrays for the dynamic programming
V = zeros(length(xc),length(tvec));
U = V; %Here big V is the value function,
%The last value of the vectors should be zero because there is no value
%after time T

small = 1e-3; %small number for so that there is no overflow in the simulation (no infinty numbers)

%new loop for the discounting, the loop changes the discounting rate for
%the single species model.
for k=1:length(dr)
    %Dynamic programming
    for i=length(tvec)-1:-1:1
        
        dVdx=ddxl*V(:,i+1);
        
        for j=1:length(dVdx)
            
            ustar(j)=1/4/max(dVdx(j)^2*xc(j),small);
            ustar(j)=max(0,min(ustar(j),Umax));
            ustar(1)=0;
        end
        
        U(:,i)=ustar;
        
        %Closed loop generator
        Gcl=G-(ustar.*xc)'.*ddxl;
        
        %euler stepping of value function
        LV=Gcl*V(:,i+1);
        V(:,i)=V(:,i+1)+LV*dt+sqrt(ustar'.*xc')*dt-dr(k)*V(:,i+1)*dt;
    end
    U0(:,k)=U(:,1);
    V0=V(:,1);
    V0=V0-min(V0);
    
    %plotting the value function and policy for each discounting rate
    figure(1)
    hold on
    plot(xc,V0,'LineWidth',2)
    xlabel('Biomass x')
    ylabel('Value V_0(x)')
    
    
    figure(2)
    hold on
    plot(xc,U0(:,k).*xc','LineWidth',2)
    xlabel('Biomass x')
    ylabel('Optimal Policy µ^*(x)')
    
    %Figure(3) plots the zoomed in poilicy plot
    figure(3)
    hold on
    plot(xc,U0(:,k).*xc','LineWidth',2)
    axis([0 1 0 2])
    xlabel('Biomass x')
    ylabel('Optimal Policy µ^*(x)')
    
end

%Labels are outside the loop
figure(1)
legend('0%','2%', '5%', '10%', '15%')
figure(2)
b = [0 0 1 1 0]; c = [0 2 2 0 0]; plot(b,c,'--','color',[0.88,0.14,0.14], 'LineWidth',2 )
legend('0%','2%', '5%', '10%', '15%')
figure(3)
b = [0 0 1 1 0]; c = [0 2 2 0 0]; plot(b,c,'--','color',[0.88,0.14,0.14], 'LineWidth',2 )
legend('0%','2%', '5%', '10%', '15%')


%Population dynamics simulations
dt=0.0001;
T=100; %time span
tvec=0:dt:T;

N=length(tvec); %intervals (how many)
Ns=1; %number of realizations, the same noise is used for every simulations



[W,Tw,dW]=ScalarStdWienerProcess(T,N,Ns); %brownian motion

%simulating the population dynamics for each discounting policy
X = zeros(size(W)); %vector for the simulation
for f=1:length(dr)
    for i=1:Ns
        X(i,1) = x0;
        for k=1:N
            dt = Tw(k+1)-Tw(k); %interval (increase)
            a=ceil(X(i,k)/0.05); %Find the nearest index in the ustar solution (how much to harvest)
            X(i,k+1) = X(i,k)+(r*X(i,k)*(1-X(i,k))-U0(a,f)*X(i,k))*dt+sigma*X(i,k)*dW(i,k);
            
        end
    end
    %plot the population dynamics
    X(end)=[];
    figure(4)
    hold on
    plot(tvec,X)
    xlabel('Time [t]')
    ylabel('Biomass')
end
figure(4)
legend('0%','2%', '5%', '10%', '15%')


