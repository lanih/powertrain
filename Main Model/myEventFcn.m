function [value, isterminal, direction] = myEventFcn(t, y, yp)
    global PTS
    
    value = (y(PTS*2+1) - y(PTS*2+2)) - 3.0;   % zero when y(5) == 0.2
    isterminal = 1;       % stop the integration
    direction = -1;       % only trigger when crossing from above (i.e. decreasing)
end