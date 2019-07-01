%Finding the MSY and the optimal profit state for the logisitc growth equation
clear all
clc
%grid and inital conditions
X=0:0.01:1;
K=1;
r=1;

%The logisitic growth equation
N=@(x) r*x.*(1-x./K);

%The surplus extracted at every state and the profit calculated
N_dt=N(X);
profit=sqrt(N(X).*X);

%The surplus and profit are then plotted
hold on
plot(X,N_dt,'LineWidth',2)
a=find(profit==max(profit));
plot(X(a),N(X(a)),'g*','LineWidth',2)

%Plotting some refrence lines, such as the state where the MSY is obtained
%and how much it produces. Inserted manually
axis([0 1 0 0.3])
line([0 0.5],[0.25 0.25],'LineWidth',1, 'Color','red','LineStyle','-.')
line([0.5 0.5],[0 0.25],'LineWidth',1, 'Color','red','LineStyle','-.')
xlabel('Biomass (x)')
ylabel('dX/dt')
