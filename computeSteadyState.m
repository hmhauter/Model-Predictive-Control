function [xs, hs] = computeSteadyState(p, F)
    xs0 = 500 * ones(4,1);
    % Initial guess of steady state
    us = [F(1); F(2)];
    % Initial Guess of steady state disturbances 
    ds = [F(3); F(4)];
    
    % Steady-state inputs
    xs = fsolve(@FourTankSystemModifiedFsolve, xs0, [], us, ds, p);
    hs = xs./ (p(12)*p(5:8));

end

