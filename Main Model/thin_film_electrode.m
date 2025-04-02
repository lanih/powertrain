function voltage = thin_film_electrode(delta)

cell=struct()

% p denotes positive (cathode), n denotes negative (anode) 
N = 10; % # of INTERNAL nodes 
h = 1/(N+1); % spacing 
Rp = 1e-6; %radius of sphere (um)
ep = 0.4; %porosity: space between particles and cell (non active material) (%) 
lp = 4.3e-5; %length of cathode 
ctp = 45829; % max concentration of li the cathode can hold 
Dp = 7.5e-14; %diffusivity term 
F = 96487; % Faraday constant 
R = 8.314; % gas constant 
soc = 0.424;%SOC: state of charge - why are these different 

iapp = delta * N*F*Dp*ctp/Rp 

%Initial condition for concentrations
%y = ones(Nx,1);
for i =1:N+2
    y(i) = soc;
end

y(N+2) = 600; % scaled - lani 

yp = zeros(N+2,1);

%time integration 
tspan = [0 3600];

%solution 
%tic;
opts = odeset('Stats','on','RelTol',1e-2,'AbsTol',1e-2,MaxOrder=5,BDF='on');
[t,y] = ode15i(@(t,y,yp) dyneqn(t,y,yp,N, iapp,Rp, Dp), tspan, y,yp,opts);
%toc;
 
 % defining the surface concentration and Voltages
 %Cathode
 ysurf = y(:,N+2); % voltage across cell at t_final 

 %Positive elctrode voltage
 %Vp = y(:,2*(Nx+2)+1);

 voltage = ysurf(N+2); % voltage at t_final at surface 
end
