function Full_SPM_Matlab_no_OCV_Tushar

global Nx delta h Rpp Rpn ep en lp ln ctp ctn Dp Dn F R socn socp ap an deltap deltan

params = struct()

params.Nx = 2;
params.h = 1/(Nx+1);
params.delta = -0.1;
Rpp = 1e-6;
Rpn = 1e-6;
ep = 0.4;
en = 0.38;
lp = 4.3e-5;
ln = 4.65e-5;
ctp = 45829;
ctn = 30555;
Dp = 7.5e-14;
Dn = 3.2e-14;
F = 96487;
R = 8.314;
socn = 0.9867;
socp = 0.424;
ap = 3*(1-ep)/Rpp;
an = 3*(1-en)/Rpn;


deltap=(-18.3*Rpp/ctp/Dp/ap/lp/F);
deltan= (18.3*Rpn/ctn/Dn/an/ln/F);

%Initial condition for concentrations
%y = ones(Nx,1);
for i =1:Nx+2
    y(i) = socp;
end

for i =Nx+3:2*(Nx+2)
    y(i) = socn;
end

yp = zeros((2*(Nx+2)),1);

%time integration 
tspan = [0 3600];

%solution 
tic;
opts = odeset('Stats','on','RelTol',1e-2,'AbsTol',1e-2,MaxOrder=5,BDF='on');
[t,y] = ode15i(@(t,y,yp) dyneqn(t,y,yp), tspan, y,yp,opts);
toc;

 % Creating the x-space
 x1_N = linspace(h,(Nx+1)*h,Nx);
 
 % defining the surface concentration
 %Cathode
 ysurfp = y(:,Nx+2);
 %Anode
 ysurfn = y(:,2*(Nx+2));
 
 figure(1);
 plot(t,[ysurfp]);
 hold on
 plot(t,[ysurfn]);
 xlabel("Time (s)");
 ylabel("Surface Concentration");