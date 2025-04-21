function [speed , rpm , motor_torque , power, voltage, current] = state_space(param, cell)

%% Net force and acceleration calculation 

F_a = param.max_torque * param.pedal_percentage / param.tire_radius; % N

net_force = F_a - param.F_d - param.F_inertia;
acceleration = net_force / param.mass; 

%% Torque is function of RPM, Voltage 

power = zeros(1,param.time_steps);
motor_torque = zeros(1,param.time_steps);
speed = zeros(1,param.time_steps);
rpm = zeros(1,param.time_steps);
voltage = zeros(1,param.time_steps);
current = zeros(1,param.time_steps); 

%% Initialize first values 

speed(1) = param.initial_speed; 
rpm(1) = param.initial_speed/param.tire_radius; 
current(1) = param.initial_motor_torque*param.initial_speed/param.tire_radius;; 
voltage(1) = param.initial_voltage;
power(1) = param.initial_motor_torque*rpm(1);
motor_torque(1) = param.initial_motor_torque;

t0 = 0;

% initial unscaled SOC from safari paper, SPM model
y0_guess  = [51555*0.5*ones(cell.battery_pts, 1); 30555*0.74*ones(cell.battery_pts,1); posocp(0.5); negocp(0.74)];   
yp0_guess = zeros(2*cell.battery_pts+2, 1);        

y_fixed   = [ones(2*cell.battery_pts, 1); zeros(2,1)];   
yp_fixed  = [zeros(2*cell.battery_pts, 1); zeros(2,1)];  


% Compute consistent initial conditions
i_guess = 10;
[y0, yp0] = decic(@(t,y,yp)spm_resf(t, y, yp, cell.battery_pts, i_guess), t0, y0_guess, y_fixed, yp0_guess, yp_fixed);

y00 = y0;
yp00 = yp0;

for t = 2:param.time_steps

    new_speed = speed(t-1) + acceleration; 
    new_rpm = new_speed/param.tire_radius;

    new_motor_torque = torque_speed_curve(new_rpm);

    power(t) = (new_motor_torque * new_rpm); % W 
    motor_torque(t) = new_motor_torque; 
    speed(t) = new_speed; 
    rpm(t) = new_rpm;

    % TODO: need to scale current
    [y0, yp0, volt] = spm_sim_statespace(t-1, y00, yp00, cell.battery_pts, current(t-1), 1.0, cell);
    voltage(t) = volt; 

    y00 = y0;
    yp00 = yp0;

    current(t) = power(t)/voltage(t); 

    
end

end



