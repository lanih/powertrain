+ clear;

PTS = 50;

% Suppose t0 = 0;
t0 = 0;
% Initial guesses
y0_guess  = [51555*0.5*ones(PTS, 1); 30555*0.74*ones(PTS,1)];   % all u[i] start at 0.1
yp0_guess = zeros(2*PTS, 1);        % assume zero initial time derivative

% Fix all values in y0_guess (1 for fixed, 0 for free)
y_fixed   = ones(2*PTS, 1);         % fix all u[i] values
yp_fixed  = zeros(2*PTS, 1);        % let ode15i compute consistent derivatives

% Compute consistent initial conditions
[y0, yp0] = decic(@battery, t0, y0_guess, y_fixed, yp0_guess, yp_fixed);

% Time span
tspan = [0 3500];

% Options with custom tolerances
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);

% Solve with ode15i
[t, y] = ode15i(@battery, tspan, y0, yp0, options);

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
