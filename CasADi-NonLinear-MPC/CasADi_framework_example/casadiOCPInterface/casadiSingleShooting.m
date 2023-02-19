function [opti, x, u] = casadiSingleShooting(model, U0, Ts, x0, D, p)
% Applies multiple shooting to solve OCP

%% Dimensions
nx = size(x0, 1); nu = size(U0, 1); nd = size(D , 1);
N  = size(U0, 2);

%% Casadi integrator
F = casadiCreateIntegrator( model, Ts, nx, nu, nd, p );

% 
sim = F.mapaccum(N);

%% Casadi optimization
% Casadi optimization object
opti = casadi.Opti();

% Variables
u  = opti.variable (nu, N)  ;       % Manipulated variables
d  = opti.parameter(nd, N)  ;       % Know disturbances

% Symbolic integration 
x = sim(x0, u, d);          
x = full([x0 x]);

% Initial condition
opti.set_initial( u, U0 );

% Casadi options - tolerance does not work as intended!
p_opts = struct;
s_opts = struct;
s_opts.acceptable_tol = 1.0e-3;
s_opts.max_iter = 5000;

% Casadi solver
opti.solver( 'ipopt', p_opts, s_opts );

% Specific values for d, p, x0
opti.set_value( d , D );

end