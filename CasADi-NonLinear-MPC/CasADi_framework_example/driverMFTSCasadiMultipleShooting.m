close all
clear
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath(genpath( 'casadiOCPInterface' ));
addpath(genpath( 'models'             ));
addpath(         'util'                );


% Set casadiPath
addpath('/home/hhauter/Documents/W22/MPC/casadi-linux-matlabR2014b-v3.5.5')
import casadi.*

%% Save Plot/Data

% Save booleans
savePlot = 0;

% Names and paths
name      = 'MFTSMultipleShooting';
figPath   = './figures/';

%% Simulation Info

% Simulation period and sampling time
t0    = 0;                      % [s]       :   initial time
tf    = 30*60;                  % [s]       :   final time
Ts    = 15.0;                   % [s]       :   sampling time
N     = (tf-t0)/Ts;             % []        :   sampling intervals
M     = 10;                     % []        :   Rk4 steps in OCP

% Initial condition
x0 = [ 5000.0 ; 5000.0 ; 5000.0; 5000.0 ];
u0 = [  200.0 ;  250.0 ];
d0 = [   40.0 ;   60.0 ];

% Boundaries
uMin = [   0.0;   0.0 ];      	% [cm^3/s]	:   input lower bound
uMax = [ 500.0; 500.0 ];        % [cm^3/s] 	:   input upper bound
dMin = [   0.0;   0.0 ];        % [cm^3/s]	:   disturbance lower bound
dMax = [ 100.0; 100.0 ];        % [cm^3/s]	:   disturbance lower bound
zMin = [20.0; 20.0; 20.0; 20.0];
%% Model

theta     = parametersMFTS();
driftFun  = @driftMFTS;
outputFun = @outputMFTS;

nx = theta.nx;
nu = theta.nu;
nd = theta.nd;
nz = theta.nz;

%% Input and Disturbance Initialisation

% Initialize
U = repmat( u0, 1, N );
D = repmat( d0, 1, N );

%% Setpoint

% Define setpoints
zBars = ...
    [ 40.0, 20.0 ;
      20.0, 25.0 ];
zBar = zeros( nz, N );
zBar(:,1:round((N)/2))     = repmat( zBars(:,1), 1, round((N)/2)   );
zBar(:,round((N)/2)+1:end) = repmat( zBars(:,2), 1, (N)-round((N)/2) );

%% Solve OCP

% Define autonomous function definition
model = @(x, u, d, p) driftFun(0, x, u, d, theta); 

% Initial guess for optimization
X0 = repmat( x0, 1, N+1 );
U0 = repmat( u0, 1, N   );
Ct = randn(1, N);
% Setup multiple shooting constraints
[opti, x, u] = casadiMultipleShooting(model, X0, U0, Ts, M, x0, D, theta);

% Output function
z = outputFun(0, x, theta);

% Objectives
alpha = 0.995;
phiz  = sum( ( reshape( z(:,2:end)  , nz*N    , 1 ) - ...
               reshape( zBar        , nz*N    , 1 ) ).^2 );
phidu = sum( ( reshape( u(:,2:end)  , nu*(N-1), 1 ) - ...
               reshape( u(:,1:end-1), nu*(N-1), 1 ) ).^2 );

% Define optimal control problem (OCP)
% ...
% objective
% opti.minimize( alpha*phiz + ( 1 - alpha )*phidu );
% constraints

%% OPTIMAL CONTROL PROBLEM 

% sla = zeros(2, N);
% objective
for i = 1:length(Ct)
    % slack variables
    s_1 = max(0, zBar(1, i) - z(1, i));
    s_1 = max(s_1, z(1,i) - zBar(1,i));


    s_2 = max(0, zBar(2, i) - z(2, i));
    s_2 = max(s_2, z(2,i) - zBar(2,i));
%     sla(:, i) = [s_1; s_2];
    opti.subject_to( 15.0 - s_1    <= z(1, i)              );
    opti.subject_to( 20.0 - s_2    <= z(2, i)              );
    opti.minimize( (Ct(i) * u(1, i) + Ct(i) * u(2, i)) + (s_1^2 + s_2^2) );
end
% constraints
opti.subject_to( uMin(1) <= u(1,:) <= uMax(1) );
opti.subject_to( uMin(2) <= u(2,:) <= uMax(2) );
opti.subject_to( 0.0     <= x(:)              );
opti.subject_to( 10.0    <= z(:)              );

% opti.subject_to(z(:, 2:end)==zBar);

% Solve OCP
sol = opti.solve();

% Extract solution
U = sol.value(u);
X = sol.value(x);
T = t0:Ts:tf;

% Simulate with optimal inputs
% [T, X] = openLoopSim( driftFun, t0:Ts:tf, x0, U, D, theta );
% X = X';

% Compute outputs
Z = outputFun(T, X, theta);

%% Plot

% Fontsize and linewidth
fs = 24;
lw = 4;

% Limits
tLim = [ T(1)   , T(end)               ]/60;
xLim = [ 0.0e+0 , 1.50*max( max( X ) ) ]/1000;
uLim = [ uMin(1), 1.25*uMax(1)         ];
dLim = [ dMin(1), 1.25*dMax(1)         ];
zLim = [ 0.0e+0 , 1.50*max( max( Z ) ) ];

% Figure
fig = figure('units','normalized','outerposition',[0.5 0 0.5 1.0]);
% tiledlayout( 2, 2, 'Padding', 'tight' );


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
stairs( T/60, [ zBar(1,:), nan ], 'linewidth', ceil(lw/2), ...
    'displayname', '$\bar{z}_{1}$', ...
    'color', 'k', 'linestyle', ':'  )
stairs( T/60, [ zBar(2,:), nan ], 'linewidth', ceil(lw/2), ...
    'displayname', '$\bar{z}_{2}$', ...
    'color', 'k', 'linestyle', '--' )

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
stairs( (t0:Ts:tf)/60, [U(1,:) U(1,end)], 'linewidth', lw, ...
    'displayname', '$F_{\mathrm{in},1}$' )
hold on
stairs( (t0:Ts:tf)/60, [U(2,:) U(2,end)], 'linewidth', lw, ...
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
stairs( (t0:Ts:tf)/60, [D(1,:) D(1,end)], 'linewidth', lw, ...
    'displayname', '$F_{d,1}$' )
hold on
stairs( (t0:Ts:tf)/60, [D(2,:) D(2,end)], 'linewidth', lw, ...
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

