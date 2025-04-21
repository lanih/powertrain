%% Car Parameters and Initial States 

param = struct(); 
cell = struct();

param.pedal_percentage = 1.0; % relative to pedal % of degree travel 

param.max_torque = 10; % Nm 
param.initial_motor_torque = 0; % Nm
param.initial_speed = 0; % m/s 
param.tire_radius = 0.4064; % m
param.F_d = 2; % N 
param.F_inertia = 4; % N
param.mass = 180; % kg
param.time_steps = 100; % s
param.initial_voltage = 4.214; % V


cell.battery_pts = 10;
% cell.radius = .08; % m^2
% cell.length = .02; % m 
cell.area = .09; 

%Cathode
cell.Rp = 2.0e-6;    
cell.Dsp = 1.0e-14;
cell.eps_p = 0.385;
cell.Lp = 80e-6;

%Anode
cell.Rn = 2.0e-6;    
cell.Dsn = 2.0e-14;
cell.eps_n = 0.485;
cell.Ln = 88e-6;

%% Call state_space function 


[speed , rpm , motor_torque , power, voltage, current] = state_space(param, cell);


%% Plots
figure(1);
subplot(3,1,1);
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

figure(3)
plot(1:param.time_steps,current*(cell.area))
xlabel( 'Time Step');
ylabel('Current (A)');
title("Current vs Time")