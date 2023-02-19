function [ h ] = outputMFTS( t, x, theta )
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
%       [ h ] = outputMFTS( t, x, theta )
%
%   Description:
%       output function for a modified four-tank system (MFTS).
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



end