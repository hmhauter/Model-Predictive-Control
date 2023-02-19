function [ xkp1k, Pkp1k ] = cdEKFPredictionERK1( driftFun, diffusionFun, ...
                        tspan, dt, xkk, Pkk, uk, dk, parametersModel )
%--------------------------------------------------------------------------
%   Author(s):
%       Marcus Krogh Nielsen
%       Addapted by: Morten Ryberg Wahlgreen
%
%   Email:
%       mkrni@dtu.dk, morwa@dtu.dk
%
% Prediction function for continuous-discrete extended Kalman filter using 
% an explicit Euler scheme for integration. 
%
%--------------------------------------------------------------------------
%% Unpack Parameters

% Number of steps
tk   = tspan(1);
tkp1 = tspan(end);
N    = (tkp1-tk)/dt;
% ...

%% State Prediction

% Time
tn = tk;
% State and covariance
xn = xkk;
Pn = Pkk;
% ...

% Main loop
for i=1:N
    
    % Compute drift term
    [ fn, An ] = driftFun(     tn, xn, uk, dk, parametersModel );
    % compute diffusion term
    [ Gn ]     = diffusionFun( tn, xn, uk, dk, parametersModel );
    % ...
   
    % Compute step
    tn = tn + dt;
    xn = xn + fn*dt;
    Pn = Pn + ( An*Pn + Pn*An' + Gn*Gn' )*dt;
    % ...
    
end
% ...

% Define one-step prediction
xkp1k = xn;
Pkp1k = Pn;
% ...

end