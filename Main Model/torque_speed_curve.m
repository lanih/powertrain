function torque = torque_speed_curve(rpm) 

if rpm < 1000
        torque = .005 * rpm + 10 ;
    elseif rpm < 5500 
        torque = 15;
    else
        torque = 15;
end

