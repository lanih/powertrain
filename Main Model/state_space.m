function [speed , rpm , motor_torque , power, voltage] = state_space(param)

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
current(1) = param.initial_current; 
voltage(1) = param.initial_voltage;

for t = 2:param.time_steps

    new_speed = speed(t-1) + acceleration; 
    new_rpm = new_speed/param.tire_radius;

    new_motor_torque = torque_speed_curve(new_rpm);

    power(t) = (new_motor_torque * new_rpm)/1000; % kW 
    motor_torque(t) = new_motor_torque; 
    speed(t) = new_speed; 
    rpm(t) = new_rpm;

    voltage(t) = thin_film_electrode(current(t-1)); 

    current(t) = power(t)/voltage(t); 

    
end




