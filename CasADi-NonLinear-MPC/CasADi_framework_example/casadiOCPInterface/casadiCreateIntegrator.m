function F = casadiCreateIntegrator( model, Ts, M, nx, nu, nd, p )
% Creates symbolic integration function that integrates Ts ahead in time

%% Casadi setup
import casadi.*
% Symbolic variables
x0 = MX.sym('x', nx, 1); u  = MX.sym('u', nu, 1); d  = MX.sym('d', nd, 1);
x  = x0;

% Symbolic right hand side (call to previous implemented function)
rhs = model(x, u, d, p);

% RK4 implementation
dt = Ts/M;
f  = Function('f', {x, u, d}, {rhs}, {'x', 'u', 'd'}, {'rhs'});
for j=1:M
    k1 = f(x            , u, d);
    k2 = f(x + dt/2 * k1, u, d);
    k3 = f(x + dt/2 * k2, u, d);
    k4 = f(x + dt   * k3, u, d);
    x  = x + dt/6*(k1 + 2*k2 +2*k3 + k4);
end
xf = x;

% Create integrator function that integrators from {x,u} to xf
F = Function('F', {x0, u, d}, {xf}, {'x', 'u', 'd'}, {'xf'});

end