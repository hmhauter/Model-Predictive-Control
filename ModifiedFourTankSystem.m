function xdot = ModifiedFourTankSystem(t,x,u,d,p)
    % MODIFIEDFOURTANKSYSTEM Model dx/dt = f(t,x,u,d,p) for Modified 4-tank System
    %
    % This function implements a differential equation model for the
    % modified 4-tank system.
    %
    % Syntax: xdot = ModifiedFourTankSystem(t,x,u,d,p)
    % Unpack states, MVs, and parameters
    m = x; % Mass of liquid in each tank [g]
    F = [u; d]; % Flow rates [cm3/s]
    a = p(1:4,1); % Pipe cross sectional areas [cm2]
    A = p(5:8,1); % Tank cross sectional areas [cm2]
    gamma = p(9:10,1); % Valve positions [-]
    g = p(11,1); % Acceleration of gravity [cm/s2]
    rho = p(12,1); % Density of water [g/cm3]
    % Inflows
    qin = zeros(4,1);
    qin(1,1) = gamma(1)*F(1); % Inflow from valve 1 to tank 1 [cm3/s]
    qin(2,1) = gamma(2)*F(2); % Inflow from valve 2 to tank 2 [cm3/s]
    qin(3,1) = (1-gamma(2))*F(2); % Inflow from valve 2 to tank 3 [cm3/s]
    qin(4,1) = (1-gamma(1))*F(1); % Inflow from valve 1 to tank 4 [cm3/s]
    % Outflows
    h = m./(rho*A); % Liquid level in each tank [cm]
    qout = a.*sqrt(2*g*h); % Outflow from each tank [cm3/s]
    % Differential equations
    xdot = zeros(4,1);
    xdot(1,1) = rho*(qin(1,1)+qout(3,1)-qout(1,1)); % Mass balance Tank 1
    xdot(2,1) = rho*(qin(2,1)+qout(4,1)-qout(2,1)); % Mass balance Tank 2
    xdot(3,1) = rho*(qin(3,1)-qout(3,1)+F(3)); % Mass balance Tank 3
    xdot(4,1) = rho*(qin(4,1)-qout(4,1)+F(4)); % Mass balance Tank 4
end