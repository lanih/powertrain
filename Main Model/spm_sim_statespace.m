function [y00,yp00,voltage] = spm_sim_statespace(t0, y0, yp0, PTS, iapp, dt)
    sol = ode15i(@(t,y,yp)spm_resf(t, y, yp, PTS, iapp), [t0 t0+dt], y0, yp0);
    t = sol.x;
    y = sol.y';

    [y00,yp00] = deval(sol,t(end));
    voltage = y(end,2*PTS+1)-y(end,2*PTS+2);
    
    %stopp = sol.xe;
    
    % if isempty(stopp) % no stop condition triggered
    %     % compute next states
    %     [y00,yp00] = deval(sol,t(end));
    % else
    %     break;
end

% function [value, isterminal, direction] = STOP_undervolt(t, y, yp, PTS)
%     value = (y(PTS*2+1) - y(PTS*2+2)) - 3.0;  
%     isterminal = 1;
%     direction = -1; % crossing from above
% end

% Options with custom tolerances
% options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MaxStep', 1, 'Events', @(t,y,yp)STOP_undervolt(t,y,yp,PTS),MaxOrder=2,BDF='on');