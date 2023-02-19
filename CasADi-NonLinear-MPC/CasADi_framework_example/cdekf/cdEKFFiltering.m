function [ xkk, Pkk, ykkm1 ] = cdEKFFiltering( tk, yk, xkkm1, Pkkm1, ...
                                        measurementFun, parametersModel )
%--------------------------------------------------------------------------
%   Author(s):
%       Marcus Krogh Nielsen
%       Addapted by: Morten Ryberg Wahlgreen
%
%   Email:
%       mkrni@dtu.dk, morwa@dtu.dk
%
% Filtering function for continuous-discrete extended Kalman filter
%
%--------------------------------------------------------------------------
%% Unpack Model Parameters

% System size
nx = length(xkkm1);
% ...

% Measurement covariance
R  = parametersModel.R;
% ...

% Identity
I = eye(nx);
% ...

%% Measurement Prediction

% Compute measurement and derivative
[ ykkm1, Ck ] = measurementFun( tk, xkkm1, parametersModel );
% ...

%% Filtering

% Innovation
ek  = yk - ykkm1;
% innovation covariance
Rek = Ck*Pkkm1*Ck' + R;
% ...

% Kalman gain
Kfxk = Pkkm1*Ck'/Rek;
% ...

% Filtered state
xkk = xkkm1 + Kfxk*ek;
% Filtered state covariance
T   = I - Kfxk*Ck;
Pkk = T*Pkkm1*T' + Kfxk*R*Kfxk';
% ...

end