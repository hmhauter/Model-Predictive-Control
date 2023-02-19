function [opti, x, u, d, x0] = casadiCollocation(model, X0, U0, Ts, M, x0, D, p)
% Applies collocation to solve OCP

%% Dimensions
nx = size(X0, 1); nu = size(U0, 1)      ; nd = size(D , 1);
N  = size(D , 2);

if (size(X0, 2)-1)/N ~= M
    error('X0 and M does not match');
end

%% Save x0 as different name (overwritten later)
x00 = x0;

%% Rename model 
f = model; 

%% Size of time intervals
dt = Ts/M;

%% Casadi optimization
% Casadi optimization object
opti = casadi.Opti();

% Variables
x  = opti.variable (nx, M*N+1);     % States
u  = opti.variable (nu, N)    ;     % Manipulated variables
d  = opti.parameter(nd, N)    ;     % Know disturbances
x0 = opti.parameter(nx, 1)    ;     % Initial state

% Constraints
opti.subject_to( x(:,1)  == x0 );      % Fix x0 

% Collocation continuity constraints (explicit euler steps)
k = 0;
for j = 1:N
    for i = 1:M
        k = k + 1;
        opti.subject_to( x(:,k+1) == ...
            x(:,k) + f(x(:,k+1), u(:,j), d(:,j), p)*dt );
    end 
end

% Initial condition
opti.set_initial( x, X0 );
opti.set_initial( u, U0 );

% options
p_opts = struct;
s_opts = struct;
% acceptable tolerance
s_opts.acceptable_tol             = 1.0e-7;
s_opts.acceptable_dual_inf_tol    = 1.0e-7;
s_opts.acceptable_constr_viol_tol = 1.0e-7;
s_opts.acceptable_compl_inf_tol   = 1.0e-7;
s_opts.acceptable_obj_change_tol  = 1.0e-7;
s_opts.acceptable_iter            = 10;
% convergence tolerance
s_opts.tol                        = 1.0e-7;
s_opts.dual_inf_tol               = 1.0e-7;
s_opts.constr_viol_tol            = 1.0e-7;
s_opts.max_iter                   = 5000;
% misc options
s_opts.print_level                = 5;
s_opts.hessian_approximation      = 'limited-memory';

% Casadi solver
opti.solver( 'ipopt', p_opts, s_opts );

% Specific values for d, p, x0
opti.set_value( x0, x00 );
opti.set_value( d, D );

end