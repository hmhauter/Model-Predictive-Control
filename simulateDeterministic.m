function [T, H, Qout] = simulateDeterministic(F, p)

% input F [cm3/s] Flow rate from pump i

t0 = 0.0; % [s] Initial time
tf = 20*60; % [s] Final time
m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
m40 = 0.0; % [g] Liquid mass in tank 4 at time t0
x0 = [m10; m20; m30; m40];
u = [F(1); F(2)];
% --------------------------------------------------------------

d = [F(3); F(4)];

% --------------------------------------------------------------
% Compute the solution / Simulate
% --------------------------------------------------------------

% Solve the system of differential equations
[T,X] = ode45(@ModifiedFourTankSystem,[t0 tf],x0,[], u,d,p);

% --------------------------------------------------------------
% help variables
[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';
rho = p(12);
g = p(11);
% Compute the measured variables
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = X(i,:)./(rho*A);
end
% Compute the flows out of each tank
Qout = zeros(nT,nX);
for i=1:nT
    Qout(i,:) = a.*sqrt(2*g*H(i,:));
end
% --------------------------------------------------------------
end

