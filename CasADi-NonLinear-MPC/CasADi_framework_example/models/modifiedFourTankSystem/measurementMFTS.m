function [ h, dhdx ] = measurementMFTS( t, x, theta )
%--------------------------------------------------------------------------
%   Author(s):
%       Morten Ryberg Wahlgreen
%
%   Email:
%       morwa@dtu.dk
%
%--------------------------------------------------------------------------
%   Call:
%       [ h ] = measurementMFTS( t, x, theta )
%
%   Description:
%       measurement function for a modified four-tank system (MFTS).
%   
%   Inputs:
%       t       :   time
%       x       : 	masses in tanks
%       theta   :   parameters
%
%   Outputs:
%       h       :   height of liquid columns
%
%--------------------------------------------------------------------------

%% Un-pack Parameters and Variables

% Variables
% ...
%
% ...
% states
m  = x(:,:);
% ...


% Parameters
% ...
%
% ...
% model parameters
A      = theta.A;
rho    = theta.rho;
% ...



%% Output Function 

% Compute Output
h = m(1:2,:)./( rho*A(1:2) );
% ...

%% x-derivative

if nargout > 1
    % Apply CasADI derivative function
    dhdx = full(theta.dydxFun( t, x ));
end

end