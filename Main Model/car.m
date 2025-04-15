%% Car Parameters and Initial States 

param = struct(); 

param.pedal_percentage = 1.0; % relative to pedal % of degree travel 

param.max_torque = 240; % Nm 
param.initial_motor_torque = 0; % Nm
param.initial_speed = 0; % m/s 
param.tire_radius = 0.4064; % m
param.F_d = 40; % N 
param.F_inertia = 30; % N
param.mass = 180; % kg
param.time_steps = 10; % s
param.initial_voltage = 600; % V
param.initial_current = 10; % A

param.battery_pts = 50;

%% Call state_space function 


[speed , rpm , motor_torque , power, voltage] = state_space(param);


%% Plots
figure(1);
subplot(3,1,1);1
plot(1:param.time_steps, speed, 'b');
xlabel('Time Step');
ylabel('Speed (m/s)');
title('Vehicle Speed');

subplot(3,1,2);
plot(1:param.time_steps, rpm, 'r');
xlabel('Time Step');
ylabel('RPM (rad/s)');
title('Motor RPM');

subplot(3,1,3);
plot(1:param.time_steps, power, 'g');
xlabel('Time Step');
ylabel('Power (kW)');
title('Motor Power Output');

figure(2)
plot(1:param.time_steps,voltage)
xlabel( 'Time Step');
ylabel('Voltage (V)');
title("Voltage vs Time")