close all; clear; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath(genpath( 'models' ));
addpath(         'util'    );

%% Save Plot/Data

% Save booleans
savePlot = 0;
saveData = 0;

% Names and paths
name      = 'MFTSOpenLoop';
dataPath  = './data/';
figPath   = './figures/';


%% Model

% Parameters, drift, and output functions
theta     = parametersMFTS();
driftFun  = @driftMFTS;
outputFun = @outputMFTS;

% System size
nx = theta.nx;
nu = theta.nu;
nd = theta.nd;
nz = theta.nz;

%% Simulation Info

% Simulation period and sampling time
t0    = 0;                      % [s]       :   initial time
tf    = 30*60;                  % [s]       :   final time
Ts    = 15.0;                   % [s]       :   sampling time
N     = (tf-t0)/Ts;             % []        :   sampling intervals

% Simulation scenario
TSim = ...                      % [s]       :   scenario (times)
    [   0,  10,  20,  30 ]*60; 
USim = ...                      % [L/s]     :   scenario (inputs)
    [ 100, 150, 300, 300 ;
      300, 200, 100, 100 ];

% Initial condition
x0 = [ 5000.0 ; 5000.0 ; 5000.0; 5000.0 ];

% Boundaries
uMin = [   0.0;   0.0 ];      	% [cm^3/s]	:   input lower bound
uMax = [ 500.0; 500.0 ];        % [cm^3/s] 	:   input upper bound
dMin = [   0.0;   0.0 ];        % [cm^3/s]	:   disturbance lower bound
dMax = [ 100.0; 100.0 ];        % [cm^3/s]	:   disturbance lower bound


%% Open-loop Simulation Inputs and Disturbances

% Initialize
U = USim;                       % [L/s]     :   scenario (inputs)
D = zeros( nd, length(TSim) );  % []        :   scenario (disturbances)


%% Open-loop Simulation

% Simulate
[T, X] = openLoopSim( driftFun, TSim, x0, U, D, theta );

% Output
Z = outputFun(T, X', theta);


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
tiledlayout( 2, 2, 'Padding', 'tight' );


% Heights
nexttile( 1 );
plot( T/60, Z(1,:), 'linewidth', lw, ...
    'DisplayName', '$z_{1}$', ...
    'color', [0, 0.4470, 0.7410], 'linestyle', '-'  )
hold on
plot( T/60, Z(2,:), 'linewidth', lw, ...
    'DisplayName', '$z_{2}$', ...
    'color', [ 0.8500, 0.3250, 0.0980], 'linestyle', '-' )
% options
grid on
ylabel('Heights [cm]')
xlim( tLim )
ylim( zLim )
legend( 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 2 )
set(gca, 'fontsize', fs)
hold off

% Masses
nexttile( 2 );
plot( T/60, X(:,1)/1000, 'linewidth', lw, 'displayname', '$x_{1}$' )
hold on
plot( T/60, X(:,2)/1000, 'linewidth', lw, 'displayname', '$x_{2}$' )
plot( T/60, X(:,3)/1000, 'linewidth', lw, 'displayname', '$x_{3}$' )
plot( T/60, X(:,4)/1000, 'linewidth', lw, 'displayname', '$x_{4}$' )
% options
grid on
ylabel('Masses [kg]')
xlim( tLim )
ylim( xLim )
legend( 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 2 )
set(gca, 'fontsize', fs)
hold off

% Inlet flows
nexttile( 3 );
stairs( TSim/60, U(1,:), 'linewidth', lw, ...
    'displayname', '$F_{\mathrm{in},1}$' )
hold on
stairs( TSim/60, U(2,:), 'linewidth', lw, ...
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


% Inlet flows
nexttile( 4 );
stairs( TSim/60, D(1,:), 'linewidth', lw, ...
    'displayname', '$F_{d,1}$' )
hold on
stairs( TSim/60, D(2,:), 'linewidth', lw, ...
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

% Save data
if saveData
    outStruct.Tms   = T;
    outStruct.Xms   = X;
    outStruct.Ums   = U;
    outStruct.Dms   = D;
    outStruct.Ts_ms = Ts;
    outStruct.t0_ms = t0;
    outStruct.tf_ms = tf;
    
    save( [ dataPath, name, '.mat' ], '-struct', 'outStruct' );
end

% Save plot
if savePlot
    exportgraphics( fig, [ figPath, name, '.png' ], 'resolution', 500 );
    close all
end



