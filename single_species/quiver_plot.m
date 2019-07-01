%Quiver plot for the Rosenzweig-MacArthur system
clear all
clc

%Construct the grid
[x,y] = meshgrid(0:0.08:1,0:0.08:1);

%initial conditions
Cmax=1;
K=1;
beta=3;
r=1;
dr=0.15;
epsilon=0.6;

%If the pareto optimality is to be plotted, then the price is changed to
%vectors
prHVec=1;%[1,2];%2; %
prCVec=1;%linspace(0,3,1000);

Feed=Cmax*beta.*x.*y./(beta.*x+Cmax);
N = r.*x.*(1-x./K)-Feed;
P = epsilon.*Feed-dr.*y;

%normalize the poopulation dynamics so the arrows in the quiver plot are
%equal in length
N_norm=N./sqrt(N.^2+P.^2);
P_norm=P./sqrt(N.^2+P.^2);

%plotting the quiver plot
figure
quiver(x,y,N_norm,P_norm,0.4)

hold on
%Plotting the equilbirium point manually
plot(0.111,0.395,'b*','LineWidth',2)
% plot(0.6667,0.3333,'b*','LineWidth',2)
axis([0 1 0 1])
xlabel('Biomass Prey')
ylabel('Biomass Predator')

%grid
x=0:0.01:1;
%sizecontrol
a=length(x);

Feed=@(x,y) Cmax*beta*x*y/(beta*x+Cmax);
prey=@(x,y) r*x*(1-x/K)-Feed(x,y);
pred=@(x,y) Feed(x,y)*epsilon-dr*y;

hold on

%adding the sustainable zone for the prey and predator to the quiver plot
for i=1:a
    for j=1:a
        %dN/dt
        N(i,j)=prey(x(i),x(j));
        %dNP/dt
        P(i,j)=pred(x(i),x(j));
        if max(N(i,j),0) && max(P(i,j),0) > 0
            plot(x(i),x(j),'r+')%,'LineWidth',1)
        end
        
    end
end

%Finding the optimal profit state for the deterministic model
maxProfit=zeros(length(prHVec),length(prCVec));
for g=1:length(prHVec)
    for h=1:length(prCVec)
        %find the highest profit
        profit=zeros(a,a);
        for i=1:a
            for j=1:a
                if max(N(i,j),0) && max(P(i,j),0) > 0
                    profit(i,j)=prHVec(g)*sqrt(N(i,j)*i)+prCVec(h)*sqrt(P(i,j)*j);
                end
            end
        end
        maxProfit(g,h)=max(max(profit));
        search=find(profit==max(max(profit)));
        [x1,y1]=ind2sub(size(profit),search);
        plot(x(x1),x(y1),'g*','LineWidth',2)
        
        fprintf('Optimal profit state: x=%2.3f, y=%2.3f',x(x1),x(y1))
        fprintf('Prey surplus produced in the state: %2.3f',prey(x(x1),x(y1)))
        fprintf('Predator surplus produced in the state: %2.3f',pred(x(x1),x(y1)))
        fprintf('The instantaneous profit in the state: %2.3f',profit(x1,y1))
        

    end
end



%plotting the population dynamics for the deterministic system

P0=1;%number of predators
N0=1;%number of preys
tspan=linspace(0,200,10000);

%Solving the deterministic Rosenzweig-MacArthur model with a built in solver
odefun=@(tspan,x)[r*x(1)*(1-x(1)/K)-(Cmax*beta*x(1)*x(2)/(beta*x(1)+Cmax));epsilon*(Cmax*beta*x(1)*x(2)/(beta*x(1)+Cmax))-dr*x(2)];
[~,NP] = ode45(odefun,tspan,[N0,P0]);

N=NP(:,1);
P=NP(:,2);
figure
plot(N,P)
xlabel('Biomass Prey')
ylabel('Biomass Predator')

figure
plot(tspan,N)
hold on
plot(tspan,P)
xlabel('Time')
ylabel('Biomass')
legend('Prey','Predator')




