%Optimal harvesting policy for a stochastic logistic growth equation

clc; clear all; close all;

%Parameters
K=1; %Carrying Capacity
r=1; %Growth rate
sigma=1/sqrt(2); %Noise in population dynamics
x0 = 1; %Initial value
% dr=0.02;%Discount rate

dt=0.0001;%0.0001; % time steps for euler method
% T=100; %time span
T=50;
tvec=0:dt:T;


%Upperbound on fisheris mortality
Umax=5;

%Discretization of state space (finite volume, for partial differential
%equations)
XMax=4; %Biomass
dx=0.05; %looping over 0.05 state step
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
U = V; %Here big V is the value function
%The last value of the vectors should be zero because there is no value
%after time T

small = 1e-3; %small number for so that there is no overflow in the simulation (no infinty numbers)

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
    V(:,i)=V(:,i+1)+LV*dt+sqrt(ustar'.*xc')*dt;
end
U0=U(:,1);
V0=V(:,1);
V0=V0-min(V0); 

figure
plot(xc,V0,'LineWidth',2)
xlabel('Biomass x')
ylabel('Value V_0(x)')
% title('Numerical solution for the value function')


figure
plot(xc,U0.*xc','LineWidth',2)
xlabel('Biomass x')
ylabel('Optimal policy µ^*(x)')
% title('Numerical solution for the optimal harvest rate')

dt=0.0001;%0.0001; % time steps for euler method
% T=100; %time span
T=100;
tvec=0:dt:T;
%Simulations
%T=100; %Final time
N=length(tvec); %intervals (how many)
Ns=1; %number of realizations

[W,Tw,dW]=ScalarStdWienerProcess(T,N,Ns); %brownian motion

X = zeros(size(W)); %vector for the smulation
Y = zeros(size(W)); %vector for the smulation
profX = zeros(size(W)); %vector for the smulation
profY = zeros(size(W)); %vector for the smulation
for i=1:Ns
    X(i,1) = x0;
    Y(i,1) = x0;
    for k=1:N
        dt = Tw(k+1)-Tw(k); %interval (increase)
        a=ceil(X(i,k)/0.05); %Find the nearest index in the ustar solution
        X(i,k+1) = X(i,k)+(r*X(i,k)*(1-X(i,k))-ustar(a)*X(i,k))*dt+sigma*X(i,k)*dW(i,k);
        profX(i,k)=sqrt(ustar(a)*X(i,k))*dt;
        %X(i,k+1) = X(i,k)+(r*X(i,k)*(1-X(i,k))-U(a,k)*X(i,k))*dt+sigma*X(i,k)*dW(i,k);
        Y(i,k+1) = Y(i,k)+(r*Y(i,k)*(1-Y(i,k))-Y(i,k)/2)*dt+sigma*X(i,k)*dW(i,k);
        profY(i,k)=sqrt(Y(i,k)/2*Y(i,k))*dt;
        %Fishing mortality ef U*x harvest ef það er bara u
    end
end
%removing the last elements to correspond with the tvec
%Nan's are removed and 
X(end)=[];
Y(end)=[];
Y(Y<0)=0;
Y(isnan(Y))=0;

%plotting the population dynamics for the single species case, with the
%optimal policy and a constant effort policy
figure
plot(tvec,X)
hold on
plot(tvec,Y)
line([tvec(1) tvec(end)],[0.67 0.67],'LineWidth',1, 'Color','red','LineStyle','-.')
xlabel('Time t')
ylabel('Biomass')
legend('\mu*(x)', 'Constant effort x/2')


%getting rid of valeus that are divided by zero
profX(end)=[];
profY(end)=[];
for i=1:length(Y)
    if Y(i)==0
        profY(i)=0;
    end
end

%profit plotting
figure
plot(tvec,profX)
hold on
plot(tvec,profY)
xlabel('Time t')
ylabel('Profit \surd{u_tX_t}')
legend('\mu*(x)', 'Constant effort x/2')
