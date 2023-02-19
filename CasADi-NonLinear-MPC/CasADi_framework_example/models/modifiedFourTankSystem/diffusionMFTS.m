function [ sigma ] = diffusionMFTS( t, x, u, d, theta )
%--------------------------------------------------------------------------
%   Author(s):
%       Morten Ryberg Wahlgreen
%
%   Email:
%       morwa@dtu.dk
%
%--------------------------------------------------------------------------
%   Call:
%       [ f ] = diffusionMFTS( t, x, u, d, theta )
%
%   Description:
%       evalutates a diffusion term for the modified four tank system 
%       (MFTS), which models variations in the inlet flow streams.
%
%   Inputs:
%       t       :     time
%       x       :     masses of tank contents
%       u       :     manipulated liquid in-flows
%       d       :     disturbance in-flows
%       theta   :     parameters
%
%   Outputs:
%       sigma   :     diffusion term
%
%--------------------------------------------------------------------------

%% Un-pack Parameters and Variables

% Variables
% ...
%
% ...
% inputs
F  = u(1:2);
% ...

% Parameters
% ...
%
% ...
% model parameters
gamma  = theta.gamma;
nLevel = theta.nLevel;
% ...

%% Diffusion Function

% Multiplicative noise on inlet flows
sigma = nLevel*[    gamma(1)*F(1)        , 0.0	                    ;   
                    0.0                  , gamma(2)*F(2)            ;    
                    0.0                  , ( 1 - gamma(2) )*F(2)    ;  
                    ( 1 - gamma(1) )*F(1), 0.0                      ];     
% ...

end