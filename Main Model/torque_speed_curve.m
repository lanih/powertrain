function torque = torque_speed_curve(rpm) 

if rpm < 1000
        torque = .005 * rpm + 45 ;
    elseif rpm < 5500 
        torque = 52;
    else
        torque = 52;
end

