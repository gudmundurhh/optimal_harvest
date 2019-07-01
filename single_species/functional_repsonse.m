%comaparing fucntional repsonses

%grid
x=0:0.001:1.5;
%initial conditions
Cmax=1;
K=1;
beta=3;
r=1;
dr=0.15;

%Calculating the intake rate for the different functional responses
y=x;
Feed=Cmax*beta.*x./(beta.*x+Cmax);

%plots
plot(x,Feed,'LineWidth',2)
hold on
plot(x,y,'LineWidth',2)
xlabel('Biomass Prey')
ylabel('Intake rate for predator')
legend('Linear functional response ','Type II functional response') 
