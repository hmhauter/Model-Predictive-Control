function [ theta ] = parametersMFTS()
%--------------------------------------------------------------------------
%   Author(s):
%       Morten Ryberg Wahlgreen
%       Marcus Krogh Nielsen
%
%   Email:
%       morwa@dtu.dk, mkrni@dtu.dk
%
%--------------------------------------------------------------------------
%   Call:
%       [ theta ] = parametersMFTS()
%
%   Description:
%       parameter function for the modified four-tank system (MFTS).
%
%   Inputs:
%
%   Outputs:
%       theta  	:   parameters
%
%--------------------------------------------------------------------------

%% Define parameters

% Size Defintions
% ...
%
% ...
% size(s)
theta.nx = 4;
theta.nu = 2;
theta.nd = 2;
theta.nw = 2;
theta.ny = 2;
theta.nz = 2;
% ...


% Parameters for the Four Tank System
% ...
%
% ...
% cross-sectional area of drain pipes
a1 = 1.2272;
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;
theta.a = [ a1; a2; a3; a4 ];
% cross-sectional area of tanks
A1 = 380.1327;
A2 = 380.1327;
A3 = 380.1327;
A4 = 380.1327;
theta.A = [ A1; A2; A3; A4 ];
% split-values of input flows
gamma1 = 0.58;
gamma2 = 0.68;
theta.gamma = [ gamma1; gamma2 ];
% gravitational accelration
theta.g   = 981;
% density of liquid
theta.rho = 1.0;
% ...

% Measurement covariance
theta.R = 1.0e0*eye(theta.ny);

% Process noise level
theta.nLevel = 0.1;


% Get driftFun and measurementFun derivatives with CasADi
import casadi.*
tCas = MX.sym('tCas',       1 , 1);
xCas = MX.sym('xCas', theta.nx, 1);
uCas = MX.sym('uCas', theta.nu, 1);
dCas = MX.sym('dCas', theta.nd, 1);
fExp    = driftMFTS(tCas, xCas, uCas, dCas, theta );
dfdxTmp = jacobian(fExp, xCas);
dfdxCas = Function('dfdxCas', {tCas, xCas, uCas, dCas}, {dfdxTmp});

yExp    = measurementMFTS( tCas, xCas, theta );
dydxTmp = jacobian(yExp, xCas);
dydxCas = Function('dydxCas', {tCas, xCas}, {dydxTmp});

% Save CasADi functions in parameter structure
theta.dfdxFun = dfdxCas;
theta.dydxFun = dydxCas;


end