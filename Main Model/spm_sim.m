clear;

PTS = 50;

t0 = 0;

% initial unscaled SOC from safari paper, SPM model
y0_guess  = [51555*0.5*ones(PTS, 1); 30555*0.74*ones(PTS,1); posocp(0.5); negocp(0.74)];   

% assume zero initial time derivative
yp0_guess = zeros(2*PTS+2, 1);        

% Fix all values in y0_guess (1 for fixed, 0 for free)
y_fixed   = [ones(2*PTS, 1); zeros(2,1)];   
yp_fixed  = [zeros(2*PTS, 1); zeros(2,1)];  

% Compute consistent initial conditions
[y0, yp0] = decic(@(t,y,yp)spm_resf(t, y, yp, PTS), t0, y0_guess, y_fixed, yp0_guess, yp_fixed);

% Time span
tspan = [0 3200];

function [value, isterminal, direction] = STOP_undervolt(t, y, yp, PTS)
    value = (y(PTS*2+1) - y(PTS*2+2)) - 3.0;  
    isterminal = 1;
    direction = -1; % crossing from above
end

% Options with custom tolerances
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MaxStep', 1, 'Events', @(t,y,yp)STOP_undervolt(t,y,yp,PTS),MaxOrder=2,BDF='on');

% Solve with ode15i
% @spm_resf
sol = ode15i(@(t,y,yp)spm_resf(t, y, yp, PTS), tspan, y0, yp0, options);

t = sol.x;
y = sol.y';

stopp = sol.xe;

if isempty(stopp);
    disp(true)
else
    disp(false)
end

[y1,yp1] = deval(sol,sol.x(end));


% Plot the results
figure;
plot(t, y(:,PTS), 'b-', 'DisplayName', 'Cathode c(t)');
hold on;
plot(t, y(:,2*PTS), 'r-', 'DisplayName', 'Anode c(t)');
grid on;
xlabel('Time t');
ylabel('Concentration');
title('Battery Solution w/ode15i');
legend show;

figure;
plot(t, y(:,2*PTS+1)-y(:,2*PTS+2), 'g-', 'DisplayName', 'Voltage v(t)');
hold on;
grid on;
xlabel('Time t');
ylabel('Voltage');
title('Battery Solution w/ode15i');
legend show;
