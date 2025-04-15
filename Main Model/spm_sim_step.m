clear;

PTS = 50;
C_RATE = 1.0;
IAPPT = 27.263836618115 * C_RATE;

t0 = 0;

% initial unscaled SOC from safari paper, SPM model
y0_guess  = [51555*0.5*ones(PTS, 1); 30555*0.74*ones(PTS,1); posocp(0.5); negocp(0.74)];   

% assume zero initial time derivative
yp0_guess = zeros(2*PTS+2, 1);        

% Fix all values in y0_guess (1 for fixed, 0 for free)
y_fixed   = [ones(2*PTS, 1); zeros(2,1)];   
yp_fixed  = [zeros(2*PTS, 1); zeros(2,1)];  


% Compute consistent initial conditions
[y0, yp0] = decic(@(t,y,yp)spm_resf(t, y, yp, PTS, IAPPT), t0, y0_guess, y_fixed, yp0_guess, yp_fixed);

function [value, isterminal, direction] = STOP_undervolt(t, y, yp, PTS)
    value = (y(PTS*2+1) - y(PTS*2+2)) - 3.0;  
    isterminal = 1;
    direction = -1; % crossing from above
end

% Options with custom tolerances
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MaxStep', 1, 'Events', @(t,y,yp)STOP_undervolt(t,y,yp,PTS),MaxOrder=2,BDF='on');

y00 = y0;
yp00 = yp0;
tstart = t0;
dt = 100;

tfull = [];
yfull = [];

while true
    if tstart < 1000
        cur_iapp = IAPPT;
    else
        cur_iapp = IAPPT*0.5;
    end

    sol = ode15i(@(t,y,yp)spm_resf(t, y, yp, PTS, cur_iapp), [tstart tstart+dt], y00, yp00, options);
    t = sol.x;
    y = sol.y';
    
    stopp = sol.xe;
    
    % Very bad for perf. 
    % TODO: Prealloc and subslice
    tfull = [tfull; t'];
    yfull = [yfull; y];

    if isempty(stopp) % no stop condition triggered
        % compute next states
        [y00,yp00] = deval(sol,t(end));
        tstart = t(end);
    else
        break;
    end
end

% Plot the results
% figure;
plot(tfull, yfull(:,PTS), 'b-', 'DisplayName', 'Cathode c(t)');
hold on;
plot(tfull, yfull(:,2*PTS), 'r-', 'DisplayName', 'Anode c(t)');
grid on;
xlabel('Time t');
ylabel('Concentration');
title('Battery Solution w/ode15i');
legend show;

figure;
plot(tfull, yfull(:,2*PTS+1)-yfull(:,2*PTS+2), 'g-', 'DisplayName', 'Voltage v(t)', 'LineWidth', 5);
hold on;
grid on;
xlabel('Time t');
ylabel('Voltage');
title('Battery Solution w/ode15i');
legend show;
