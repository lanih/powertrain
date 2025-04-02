%1D Slab Electrode Model 

%F5 to run the model 

clear all
clc
close all
format long

% parameter set
global Nx delta h
Nx = 4;
h = 1/(Nx-2);
delta = -0.01;

% initial condition Anode    state of charge SOC as a function of Li-ion
% concentration.
y = ones(2*Nx,1);
for i =1:Nx
y(i)=0.9867;
end
for i =Nx+1:2*Nx
y(i)=0.424;
end

yp = zeros(2*Nx,1);
for i =1:Nx
yp(i)=0.9867;
end
for i =Nx+1:2*Nx
yp(i)=0.424;
end

% integration time
tspan = [0 60];

tic;
opts = odeset('Stats','on','RelTol',1e-5,'AbsTol',1e-6,MaxOrder=2,BDF='on');
[t,y] = ode15i(@(t,y,yp) dyneqn(t,y,yp), tspan, y,yp,opts);
toc;
 
ysurf1 = y(:,Nx)/2+y(:,Nx-1)/2;
ysurf2 = y(:,Nx*2)/2+y(:,(Nx*2)-1)/2;
figure(1);
plot(t,[y(:,2:Nx-1),ysurf1])
legend('C2','C3','Csurf')
xlabel('t') 
ylabel('C') 
figure(2);
hold on
% plot(x1_N,[y(end,Nx+2:2*Nx-1)])
plot(t,ysurf1)
plot(t,ysurf2)
xlabel('t') 
ylabel('C') 

function res = dyneqn(t,y,yp)
global Nx h delta 

%Anode 
 res = zeros(2*Nx,1); %Expand row to size 2Nx
        res(1)=y(1)-y(2);  %dc/dx = 0
    for i=2:Nx-1
        %Rectangular 
        %res(i) = yp(i)-(y(i+1)-2*y(i)+y(i-1))/h^2; %dc/dt -d^c/dx^2 
        %Spherical
        %res(i) = yp(i)-((y(i-1)-2*y(i)+y(i+1))/h^2)+(2/(((i-2)*h)+(h/2))*(y(i)-y(i-1))/h); 
        res(i) = yp(i)-((y(i+1)-2*y(i)+y(i-1))/((h)^2)+(2/(h*(i-1.5))*(y(i)-y(i-1))/h)) ;

    end

 res(Nx)=(y(Nx)-y(Nx-1))/h-delta; % dc/dx = -delta Anode BC


%Cathode 
         % res = zeros(Nx+1,1);
        res(Nx+1)=y(Nx+1)-y(Nx+2);  %dc/dx = 0
    for i= Nx+2:2*Nx-1 %2:Nx-1  % Add Nx to i Nx+2:2*Nx-1
        %Rectangular 
        %res(i) = yp(i)-(y(i+1)-2*y(i)+y(i-1))/h^2; %dc/dt -d^c/dx^2 
        %Spherical
        %res(i) = yp(i)-((y(i-1)-2*y(i)+y(i+1))/h^2)+(2/(((i-2)*h)+(h/2))*(y(i)-y(i-1))/h); 
        res(i) = yp(i)-((y(i+1)-2*y(i)+y(i-1))/((h)^2)+(2/(h*(i-1.5))*(y(i)-y(i-1))/h)) ;
    end

res(2*Nx)=((y(2*Nx)-y(2*Nx-1))/h)+delta; % Cathode BC
end 





