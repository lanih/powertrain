function [y00,yp00,voltage] = spm_sim_statespace(t0, y0, yp0, PTS, iapp, dt, cell)
    
    current_density = iapp / (cell.area);
    
    
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'MaxStep', 1, 'Events', @(t,y,yp)STOP_undervolt(t,y,yp,PTS),MaxOrder=2,BDF='on');
    sol = ode15i(@(t,y,yp)spm_resf(t, y, yp, PTS, current_density, cell), [t0 t0+dt], y0, yp0, options);
    t = sol.x;
    y = sol.y';

    [y00,yp00] = deval(sol,t(end));
    voltage = y(end,2*PTS+1)-y(end,2*PTS+2);
    

end
