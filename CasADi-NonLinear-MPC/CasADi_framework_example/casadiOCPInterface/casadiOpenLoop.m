function [T, X] = casadiOpenLoop( model, tspan, Ts, x0, U, D, p )
% Performs an open loop simulation

%% Dimensions
nx = size(x0, 1);
nu = size(U , 1);
nd = size(D , 1);
N  = size(U , 2);

%% Time
t0 = tspan(1);
tf = tspan(end);

%% Casadi integrator
F = casadiCreateIntegrator( model, Ts, nx, nu, nd, p );

%% Simulate

% Create recursive call to F to make simulation function
sim = F.mapaccum(N);

% Simulate
T = linspace(t0, tf, N+1)';  
X = sim(x0, U, D);          
X = full([x0 X]);

end