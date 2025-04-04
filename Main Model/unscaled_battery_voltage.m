clear;

PTS = 50;

% Suppose t0 = 0;
t0 = 0;
% Initial guesses

% all u[i] start at 0.1
y0_guess  = [51555*0.5*ones(PTS, 1); 30555*0.74*ones(PTS,1)];   
y0_guess  = [y0_guess; 4.234963004675769; 0.09280796340471076]; 

% assume zero initial time derivative
yp0_guess = zeros(2*PTS+2, 1);        

% Fix all values in y0_guess (1 for fixed, 0 for free)
y_fixed   = [ones(2*PTS, 1); zeros(2,1)];         % fix all u[i] values
yp_fixed  = [zeros(2*PTS, 1); zeros(2,1)];       % let ode15i compute consistent derivatives

% Compute consistent initial conditions
[y0, yp0] = decic(@battery_voltage, t0, y0_guess, y_fixed, yp0_guess, yp_fixed);

% Time span
tspan = [0 3500];

% Options with custom tolerances
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @myEventFcn);

% Solve with ode15i
[t, y] = ode15i(@battery_voltage, tspan, y0, yp0, options);

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
