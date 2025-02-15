%Defining function for the discretization
%ysurf(end)
function res = dyneqn(t,y,yp,u,ks)
global Nx h deltap deltan Dp Rpp ctp lp F Dn Rpn ctn ln ap an
res = zeros(2*(Nx+2),1);
%cathode equations
res(1)=y(1)-y(2);
for i=2:Nx+1
    res(i)=yp(i)-Dp/Rpp^2*(2/(i*h)*(y(i+1)-y(i))/h+(y(i+1)-2*y(i)+y(i-1))/h^2);
end
res(Nx+2)=(y(Nx+2)-y(Nx+1))/h-(18.3*Rpp/ctp/Dp/ap/lp/F);
%Anode equations
res(Nx+3)=y(Nx+3)-y(Nx+4);
for i=Nx+4:(2*(Nx+2)-1)
    res(i)=yp(i)-Dn/Rpn^2*(1/(i*h)*(y(i+1)-y(i))/h+(y(i+1)-2*y(i)+y(i-1))/h^2);
end
res(2*(Nx+2))=(y(2*(Nx+2))-y(2*(Nx+2)-1))/h-(-18.3*Rpn/ctn/Dn/an/ln/F);

end
