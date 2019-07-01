%Testing where the numerical solution is stable by observing visually where
%the upperbound of the fisheries mortality rate and upper bound of the
%biomass are changed.

clc; clear all; close all;

%Parameters
K=2; %Carrying Capacity
r=1; %Growth rate
sigma=1/sqrt(2); %Noise in population dynamics
x0 = 1; %Initial value
dr=0.05;
dt=0.0001;% time steps for euler method
T=100; %time span

tvec=0:dt:T;

% Functions entering into the model
%Uncontrolled system

f = @(x) r.*x.*(1-x./K);

D = @(x) 1/2*sigma.^2*x.^2;
dD = @(x) sigma.^2*x;

%The optimal policy and value function are tested for each of the cases
%defined in UMaxVec and XmaxVec, and compared.
UMaxVec=[3,6,10];
XmaxVec=[5,6,7];


for e=1:length(XmaxVec)
    for g=1:length(UMaxVec)
        %Upperbound on fisheris mortality
        Umax=UMaxVec(g);
        
        %Discretization of state space
        XMax=XmaxVec(e); %Biomass
        dx=0.05; %looping over 0.05 time steps
        xi=0:dx:XMax;
        xc=xi(1:length(xi)-1)+0.025;
        
        
        % Functions entering into the model, has to be updated
        %Uncontrolled system
        v = @(x) f(x)-dD(x);
        G = fvade(v,D,xi,'r');
        
        % Effect of the fishing: The "generator" d/dx %clarification again
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
            
            %euler stepping of value function)
            LV=Gcl*V(:,i+1);
            V(:,i)=V(:,i+1)+LV*dt+sqrt(ustar'.*xc')*dt-dr*V(:,i+1)*dt; 
            U0=U(:,1);
            V0=V(:,1);
            V0=V0-min(V0);
        end
        %plotting the value function and optimal policy for each case
        figure(1)
        hold on
        plot(xc,V0,'LineWidth',2)
        xlabel('Biomass x')
        ylabel('Value V_0(x)')        
        
        figure(2)
        hold on
        plot(xc,U0.*xc','LineWidth',2)
        xlabel('Biomass x')
        ylabel('Optimal Policy µ^*(x)')

    end
end

