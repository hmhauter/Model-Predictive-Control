function [ f, dfdx ] = driftMFTS( t, x, u, d, theta )
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
%       [ f ] = driftMFTS( t, x, u, d, theta )
%
%   Description:
%       evalutates the drift term of the modified four tank system (MFTS).
%
%   Inputs:
%       t       :     time
%       x       :     masses of tank contents
%       u       :     manipulated liquid in-flows
%       d       :     disturbance in-flows
%       theta   :     parameters
%
%   Outputs:
%       f       :     state drift function
%
%--------------------------------------------------------------------------

%% Un-pack Parameters and Variables

% Variables
% ...
%
% ...
% masses in containers
m  = x(1:4);
% inputs
F  = u(1:2);
% disturbances
FD = d(1:2);
% ...


% Parameters
% ...
%
% ...
% model parameters
a      = theta.a;
A      = theta.A;
gamma  = theta.gamma;
g      = theta.g;
rho    = theta.rho;
% ...



%% Drift Function

% Compute flows out of the system
% ...
%
% ...
% compute heights from masses and tank geometry
h    = m./( rho*A );
% compute out-flows
qOut = a.*sqrt( 2*g*h );
% ...


% Compute flows into the system
% ...
%
% ...
% compute in-flows
qIn = ...
    [         gamma(1)  *F(1)         + qOut(3)	;    	% flow into tank 1
              gamma(2)  *F(2)         + qOut(4)	;       % flow into tank 2
        ( 1 - gamma(2) )*F(2) + FD(1)           ;       % flow into tank 3
        ( 1 - gamma(1) )*F(1) + FD(2)           ];      % flow into tank 4
% ...


% Evaluate drift function
f = rho.*( qIn - qOut );
% ...


%% x-derivative
if nargout > 1
    % Apply CasADI derivative function
    dfdx = full(theta.dfdxFun( t, x, u, d ));
end



end