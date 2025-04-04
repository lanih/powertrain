function cell_voltage = cell_model(delta)
global Nx delta h Rpp Rpn ep en lp ln ctp ctn Dp Dn F R socn socp ap an deltap deltan

cell=struct()

% p denotes positive (cathode), n denotes negative (anode) 
Nx = 2; % # of INTERNAL nodes 
h = 1/(Nx+1); % spacing 
%delta = -1; %current density 
Rpp = 1e-6; %radius of sphere (um)
Rpn = 1e-6;
ep = 0.4; %porosity: space between particles and cell (non active material) (%) 
en = 0.38;
lp = 4.3e-5; %length of cathode 
ln = 4.65e-5;
ctp = 45829; % max concentration of li the cathode can hold 
ctn = 30555; 
Dp = 7.5e-14; %diffusivity term 
Dn = 3.2e-14;
F = 96487; % Faraday constant 
R = 8.314;
socn = 0.9867; %SOC: state of charge - why are these different 
socp = 0.424;
ap = 3*(1-ep)/Rpp; % reformulation of parameters, SA per volume ratio, want to be higher 
an = 3*(1-en)/Rpn;

F = @(t, y, yp) [
    yp(1) + y(1) - cos(t);
    yp(2) - y(1) + sin(t)
];


% deltap=delta-(-18.3*Rpp/ctp/Dn/an/ln/F); % reformulation 
% deltan= delta-(-18.3*Rpp/ctp/Dp/ap/lp/F); 

deltan=(-18.3*Rpp/ctp/Dn/an/ln/F); % reformulation 
deltap= delta * (-18.3*Rpp/ctp/Dp/ap/lp/F); 

%Initial condition for concentrations
%y = ones(Nx,1);
for i =1:Nx+2
    y(i) = socp;
end

y(Nx+1) = 1; % scaled - lani 

yp = zeros(Nx+2,1);

%time integration 
tspan = [0 3600];
iapp = delta; 
%solution 
%tic;
opts = odeset('Stats','on','RelTol',1e-2,'AbsTol',1e-2,MaxOrder=5,BDF='on');

[t,y] = ode15i(@(t,y,yp) dyneqn(t,y,yp,Nx, iapp, Rpp, Dp), tspan, y,yp,opts);
%toc;

 % Creating the x-space
 x1_N = linspace(h,(Nx+1)*h,Nx);
 
 % defining the surface concentration and Voltages
 %Cathode
 ysurfp = y(:,Nx+2);

 %Positive elctrode voltage
 %Vp = y(:,2*(Nx+2)+1);

 cell_voltage = ysurfp(Nx+2);
 
 % %Plotting surface concentrations 
 % figure(1);
 % plot(t,ysurfp);
 % hold on
 % plot(t,ysurfn);
 % xlabel("Time (s)");
 % ylabel("Surface concentration");
 % 
 % %plotting voltages
 % figure(2)
 % plot(t,Vp, "r");
 % hold on
 % plot(t,Vn, "k");
 % xlabel("Time (s)");
 % ylabel("Voltage (V)");
 % 
 % figure(3)
 % plot(t,cell_voltage);
 % xlabel("Time (s)");
 % ylabel("Voltage (V)");

end


