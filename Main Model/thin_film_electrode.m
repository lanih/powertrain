%function voltage = thin_film_electrode(delta)

delta = 200; %A/m^2
cell=struct()

N = 10; % # of INTERNAL nodes 
h = 1/(N+1); % spacing 
Rp = 1e-6; %radius of sphere (um)
ep = 0.4; %porosity: space between particles and cell (non active material) (%) 
lp = 4.3e-5; %length of cathode 
ctp = 45829; % max concentration of li the cathode can hold 
Dp = 7.5e-14; %diffusivity term 
F = 96487; % Faraday constant 
R = 8.314; % gas constant 
soc = 0.424;%SOC: state of charge  

iapp = -delta * N*F*Dp*ctp/Rp;

%Initial condition for concentrations
y_ic = ones(N+2,1);

for i =1:N+2
    y_ic(i) = soc;
end

y_ic(N+2) = 4.2; % scaled - lani 

yp_ic = zeros(N+2,1);

%time integration 
tspan = [0 3600];

%solution 

[t,y] = ode15i(@(t,y,yp) dyneqn(t,y,yp,N, iapp,Rp, Dp), tspan,y_ic,yp_ic);
 
 % defining the surface concentration and Voltages
 %Cathode
 node_voltage = y(:,N+2) % voltage across cell at t_final 

 %Positive elctrode voltage
 %Vp = y(:,2*(Nx+2)+1);

% %Plotting surface concentrations 
 figure(1)
 plot(N,node_voltage);

% figure (2)
% plot(t,voltage); 


