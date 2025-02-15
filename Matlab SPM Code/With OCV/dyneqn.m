%Defining function for the discretization
%ysurf(end)
function res = dyneqn(t,y,yp,u,ks)
global Nx h deltap deltan Dp Rpp ctp lp F Dn Rpn ctn ln ap an
res = zeros(2*(Nx+2)+2,1); 
%cathode equations
res(1)=y(2)-y(1);
for i=2:Nx+1 % Starts at 2 and goes to Nx+1 for internal nodes 
    res(i)=yp(i)-Dp/Rpp^2*(2/((i-1)*h)*(y(i)-y(i-1))/h+(y(i+1)-2*y(i)+y(i-1))/h^2);
end
res(Nx+2)=-(y(Nx+2)-y(Nx+1))/h-(deltap);
%Anode equations
res(2*Nx+1)=y(2*Nx+1)-y(2*Nx + 2); % starts at 2Nx+1 
for i=2*Nx+2:(2*Nx+3)
    res(i)=yp(i)-Dn/Rpn^2*(2/((i-1)*h)*(y(i+1)-y(i))/h+(y(i+1)-2*y(i)+y(i-1))/h^2);
end
res(2*(Nx+2))=(y(2*(Nx+2))-y(2*(Nx+2)-1))/h+(deltan);
%OCV Expression
res(2*(Nx+2)+1)=y(2*(Nx+2)+1)-(9.41894182550-20.8089294802*y(Nx+2)+26.0321098756*y(Nx+2)^2-11.2415464828*y(Nx+2)^3-292.586540767*(y(Nx+2)-.96)^2*(1/2*tanh(1000*y(Nx+2)+1000)+1/2*tanh(1000*y(Nx+2)-960.00))...
    -.212804379260*(y(Nx+2)-.733)^3*(1-1/2*tanh(1000*y(Nx+2)+1000)-1/2*tanh(1000*y(Nx+2)-733.000))+125.437252178*(y(Nx+2)-.575)^3*(1-1/2*tanh(1000*y(Nx+2)+1000)-1/2*tanh(1000*y(Nx+2)-575.000))...
    -198.413506247*(y(Nx+2)-.535)^3*(1-1/2*tanh(1000*y(Nx+2)+1000)-1/2*tanh(1000*y(Nx+2)-535.000))+164.029238273*(y(Nx+2)-.466)^3*(1-1/2*tanh(1000*y(Nx+2))-1/2*tanh(1000*y(Nx+2)-466.000)));

res(2*(Nx+2)+2)=y(2*(Nx+2)+2)-(0.7222+0.1387*y(2*(Nx+2))+0.029*y(2*(Nx+2))^(1/2)-0.0172./y(2*(Nx+2))+0.0019./y(2*(Nx+2))^1.5+0.2808*exp(0.9-15*y(2*(Nx+2)))-0.7984*exp(0.4465*y(2*(Nx+2))-0.4108));

end
