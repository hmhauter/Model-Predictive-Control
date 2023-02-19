close all
clear
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath(genpath( 'casadiOCPInterface' ));
addpath(genpath( 'models'             ));
addpath(         'util'                );
addpath(         'cdekf'               );


% Set casadiPath
addpath('/home/hhauter/Documents/W22/MPC/casadi-linux-matlabR2014b-v3.5.5')
import casadi.*

%% Save Plot/Data

% Save booleans
savePlot = 0;

% Names and paths
name      = 'MFTSCollocationClosedLoop';
figPath   = './figures/';

%% Simulation Info

% for report: allowd to use code but try diff parameter combinations etc.

% Simulation period and sampling time
t0    = 0;                      % [s]       :   initial time
tf    = 60*60;                  % [s]       :   final time
Ts    = 15.0;                   % [s]       :   sampling time
Nsim  = 10;                     % []        :   euler-maruyama steps in simulation
N     = 4.0*10.0;               % []        :   discrete control horizon 
Ns    = (tf-t0)/Ts;             % []        :   sampling intervals
M     = 10;                     % []        :   euler steps in OCP

% Initial condition
x0 = [ 5000.0 ; 5000.0 ; 5000.0; 5000.0 ];
u0 = [  200.0 ;  250.0 ];
d0 = [   40.0 ;   60.0 ];

% Boundaries
uMin = [   0.0;   0.0 ];      	% [cm^3/s]	:   input lower bound
uMax = [ 500.0; 500.0 ];        % [cm^3/s] 	:   input upper bound
dMin = [   0.0;   0.0 ];        % [cm^3/s]	:   disturbance lower bound
dMax = [ 100.0; 100.0 ];        % [cm^3/s]	:   disturbance lower bound

%% Model

% Simulation model 
pSim                = parametersMFTS();
driftFunSim         = @driftMFTS;
diffusionFunSim     = @diffusionMFTS;
measurementFunSim   = @measurementMFTS;

% Controller model (can be different from the simulation model)
pCtrl           = parametersMFTS();
driftFun        = @driftMFTS;
diffusionFun    = @diffusionMFTS;
outputFun       = @outputMFTS;
measurementFun  = @measurementMFTS;

% Dimensions (note: these can differ from simulation to controller model)
nx = pSim.nx;
nz = pSim.nz;
nd = pSim.nd;
nu = pSim.nu;
ny = pSim.ny;
nw = pSim.nw;

%% Input and Disturbance Initialisation

% Initialize
D = repmat( d0, 1, Ns+N );

%% Setpoint


% Define setpoints
zBar = zeros( nz, Ns+N*M );
zBar(1, 1:round(Ns/2) ) = 50.0;
zBar(2, 1:round(Ns/2) ) = 30.0;

zBar(1, round(Ns/2)+1:end ) = 35.0;
zBar(2, round(Ns/2)+1:end ) = 45.0;

%% Solve OCP

% Define autonomous function definition
model = @(x, u, d, p) driftFun(0, x, u, d, pCtrl); 

% Initial guess for optimization
X0 = repmat( x0, 1, M*N+1 );
U0 = repmat( u0, 1, N     );

% Setup multiple shooting constraints
[opti, x, u, d, x0Cas] = casadiCollocation(model, X0, U0, Ts, M, x0, ...
                                                        D(:,1:N), pCtrl);
% Output function
z       = outputFun(0, x, pCtrl);           % Compute symbolic output
zbar    = opti.parameter(nz, size(z,2));    % Make output target symbolic

% 
ukm1Cas = opti.parameter(nu, 1);            % Previous input

% Objectives (output target tracking and input rate of movement penalty)
alpha = 0.995;
phiz  = sum( ( reshape( z           , nz*(M*N+1), 1 ) - ...
               reshape( zbar        , nz*(M*N+1), 1 ) ).^2 );
phidu = sum((u(:, 1) - ukm1Cas).^2) + ...
            sum( ( reshape( u(:,2:end)  , nu*(N-1), 1 ) - ...
               reshape( u(:,1:end-1), nu*(N-1), 1 ) ).^2 );

% Define optimal control problem (OCP)
% ...
% objective
opti.minimize( alpha*phiz + ( 1 - alpha )*phidu );
% constraints
opti.subject_to( uMin(1) <= u(1,:) <= uMax(1) );
opti.subject_to( uMin(2) <= u(2,:) <= uMax(2) );
opti.subject_to( 0.0     <= x(:)              );

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
s_opts.hessian_approximation      = 'exact';    %'limited-memory';

% make silent
s_opts.print_level                = 0;
p_opts.print_time                 = 0;

% Casadi solver
opti.solver( 'ipopt', p_opts, s_opts );


%% Setup noise
% Seed
seed = 5;
rng(seed);

% Measurement noise
Rv = pCtrl.R;
try 
    L = chol(Rv, 'lower');
catch
    L = zeros(size(Rv));
end
v = L*randn( ny, Ns+1 );

% Process noise
dt = Ts/Nsim; % Nsim euler-maruyama in simulation (and euler steps in CDEKF)
dW = sqrt(dt) * randn(nw, Ts/dt, Ns);

%% Closed-loop

% Save for plotting
T = zeros(Ns*Nsim+1, 1);
T(1) = t0;

X = zeros( nx, Ns*Nsim+1 );
X(:, 1) = x0;

U = zeros( nu, Ns   );
Y = zeros( ny, Ns+1 );

% Simulation data
tk   = t0;
xk   = x0;

% CDEKF states
xkm1km1 = x0;
Pkm1km1 = 1.0e-5*eye(nx);

% Variables to save previous inputs and disturbances
ukm1 = u0;
dkm1 = d0;
kalman_x = [];
kalman_x_pred = [];
y_1 = [];
y_2 = [];
% Closed-loop
for k = 1:Ns
    fprintf("Iteration %4d of %4d... \n", k, Ns);
    
    % ---------------------- NMPC ---------------------- %
    % Measurement 
    yk = measurementFunSim(tk, xk, pSim) + v(:, k);

    % State estimation
    tspanCDEKF = [ tk-Ts, tk ];
    [ xkkm1, Pkkm1 ] = cdEKFPredictionERK1( driftFun, diffusionFun, ...
                    tspanCDEKF, dt, xkm1km1, Pkm1km1, ukm1, dkm1, pCtrl );
    
    [ xkk, Pkk, ykkm1 ] = cdEKFFiltering( tk, yk, xkkm1, Pkkm1, ...
                                        measurementFun, pCtrl );

    kalman_x = [kalman_x; xkk' ];
    kalman_x_pred = [kalman_x_pred; xkkm1' ];

    y_1 = [y_1; yk'];
    y_2 = [y_2; ykkm1'];

    % Save for next iteration
    xkm1km1 = xkk;
    Pkm1km1 = Pkk;

    % Disturbance in horizon
    Dk = D(:, k:k+N-1);

    % Target in horizon
    zBark = zBar(:, k:k+N*M);

    % Update CasADi OCP
    tic
    opti.set_value( x0Cas  , xkk   );       % Initial condition
    opti.set_value( d      , Dk    );       % Disturbance
    opti.set_value( zbar   , zBark );       % Target
    opti.set_value( ukm1Cas, ukm1  );       % Previous inputs

    opti.set_initial( x, X0 );          % State initial condition for OCP
    opti.set_initial( u, U0 );          % Input initial condition for OCP
    timeCasadiSetup = toc;

    % Solve OCP
    tic
    sol = opti.solve();
    timeCasadiSolve = toc;

    % Extract solution
    Xsol = sol.value(x);
    Usol = sol.value(u);

    % Implement only first inputs
    uk   = Usol(:, 1);              % Solution to implement
    ukm1 = uk;                      % Save as previous implemented inputs

    % Update initial guess
    X0 = [ Xsol(:, 2:end), Xsol(:, end) ];
    U0 = [ Usol(:, 2:end), Usol(:, end) ];

 
    % -------------------------------------------------- %
    



    % ------------------- Simulation ------------------- %
    % Current disturbances
    dk = D(:, k);

    % Simulation time span
    tspan = [tk tk+Ts];

    % Simulate system
    [Tsim, Xsim] = eulerMaruyama(driftFunSim, diffusionFunSim, ...
                            tspan, dt, xk, uk, dk, pSim, dW(:,:,k));

    % Save data
%     T = [T; Tsim];
%     X = [X, Xsim(:, 2:end)];

    T(   (2:Nsim+1) + (k-1)*Nsim ) = Tsim(2:end);
    X(:, (2:Nsim+1) + (k-1)*Nsim ) = Xsim(:,2:end);
    U(:, k) = uk;
    Y(:, k) = yk;

    % Update time and states for next iteration
    tk = Tsim(   end);
    xk = Xsim(:, end);

    % ------------------------------------------------ %

    fprintf("    CasADi setup time: %6.4f [s].\n", timeCasadiSetup);
    fprintf("    CasADi solve time: %6.4f [s].\n", timeCasadiSolve);
    fprintf("Done.\n\n");

end

% Last measurement
Y(:, end) = measurementFunSim(tk, xk, pSim) + v(:, end);


%% Plot

Z = outputFun( 0, X, pSim );

% Fontsize and linewidth
fs = 24;
lw = 6;
ms = 10;

% Limits
tLim = [ T(1)   , T(end)               ]/60;
xLim = [ 0.0e+0 , 1.50*max( max( X ) ) ]/1000;
uLim = [ uMin(1), 1.25*uMax(1)         ];
dLim = [ dMin(1), 1.25*dMax(1)         ];
zLim = [ 0.0e+0 , 1.50*max( max( Z ) ) ];

% Figure
fig = figure('units','normalized','outerposition',[0.5 0 1.0 1.0]);
% tiledlayout( 2, 2, 'Padding', 'tight' );

t = t0:Ts:tf;

% Heights
% nexttile( 1 );
subplot(2, 2, 1);
plot( T/60, Z(1,:), 'linewidth', lw, ...
    'DisplayName', '$z_{1}$', ...
    'color', [0, 0.4470, 0.7410], 'linestyle', '-'  )
hold on
plot( T/60, Z(2,:), 'linewidth', lw, ...
    'DisplayName', '$z_{2}$', ...
    'color', [ 0.8500, 0.3250, 0.0980], 'linestyle', '-' )
stairs( t/60, [ zBar(1,1:Ns), nan ], 'linewidth', ceil(lw/2), ...
    'displayname', '$\bar{z}_{1}$', ...
    'color', 'k', 'linestyle', ':'  )
stairs( t/60, [ zBar(2,1:Ns), nan ], 'linewidth', ceil(lw/2), ...
    'displayname', '$\bar{z}_{2}$', ...
    'color', 'k', 'linestyle', '--' )
plot( t/60, Y(1,:), '.', 'markersize', ms, ...
    'color', 'r', 'displayname', '$y_1$')
plot( t/60, Y(2,:), '.', 'markersize', ms, ...
    'color', 'b', 'displayname', '$y_2$')

% options
grid on
ylabel('Heights [cm]')
xlim( tLim )
ylim( zLim )
legend( 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 2 )
set(gca, 'fontsize', fs)
hold off

% Masses
% nexttile( 2 );
subplot(2, 2, 2);
plot( T/60, X(1,:)/1000, 'linewidth', lw, 'displayname', '$x_{1}$' )
hold on
plot( T/60, X(2,:)/1000, 'linewidth', lw, 'displayname', '$x_{2}$' )
plot( T/60, X(3,:)/1000, 'linewidth', lw, 'displayname', '$x_{3}$' )
plot( T/60, X(4,:)/1000, 'linewidth', lw, 'displayname', '$x_{4}$' )
% options
grid on
ylabel('Masses [kg]')
xlim( tLim )
ylim( xLim )
legend( 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 2 )
set(gca, 'fontsize', fs)
hold off

% Inlet flows
% nexttile( 3 );
subplot(2, 2, 3);
stairs( t/60, [U(1,:) U(1,end)], 'linewidth', lw, ...
    'displayname', '$F_{\mathrm{in},1}$' )
hold on
stairs( t/60, [U(2,:) U(2,end)], 'linewidth', lw, ...
    'displayname', '$F_{\mathrm{in},2}$' )
plot( [ t0, tf ]/60, [ uMin(1), uMin(1) ], 'k--', 'linewidth', ceil(lw/2), ...
    'displayname', 'bounds' )
p = plot( [ t0, tf ]/60, [ uMax(1), uMax(1) ], 'k--', ...
    'linewidth', ceil(lw/2) );
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% options
grid on
ylabel('Inlet flow rates [cm$^3$/s]')
xlabel('Time [min]')
xlim( tLim )
ylim( uLim )
legend( 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 3 )
set(gca, 'fontsize', fs)
hold off


% Disturbance flows
% nexttile( 4 );
subplot(2, 2, 4);
stairs( t/60, D(1,1:Ns+1), 'linewidth', lw, ...
    'displayname', '$F_{d,1}$' )
hold on
stairs( t/60, D(2,1:Ns+1), 'linewidth', lw, ...
    'displayname', '$F_{d,2}$' )
plot( [ t0, tf ]/60, [ dMin(1), dMin(1) ], 'k--', 'linewidth', ceil(lw/2), ...
    'displayname', 'bounds' )
p = plot( [ t0, tf ]/60, [ dMax(1), dMax(1) ], 'k--', ...
    'linewidth', ceil(lw/2) );
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% options
grid on
ylabel('Disturbance flow rates [cm$^3$/s]')
xlabel('Time [min]')
xlim( tLim )
ylim( dLim )
legend( 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 3 )
set(gca, 'fontsize', fs)
hold off


%% Save Data and Illustration

% Save plot
if savePlot
    exportgraphics( fig, [ figPath, name, '.png' ], 'resolution', 500 );
    close all
end

