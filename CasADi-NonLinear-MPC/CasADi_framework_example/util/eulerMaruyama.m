function [T, X] = eulerMaruyama(driftFun, diffusionFun, tspan, dt, x0, u, d, p, dw)
%
%   eulerMaruyama
%       Call: [T, X] = eulerMaruyama(driftFun, diffusionFun, tspan, dt, x0, u, d, p, dw)
%
%%
% Number of states
nx = length(x0);

% Number of steps
N = round((tspan(2)-tspan(1))/dt);

% Check that dw has correct size
if size(dw,2) < N
    error( 'Not enough noise realizations in dw' );
end

T = zeros(N+1, 1);
X = zeros(nx, N+1);

T(1)   = tspan(1);
X(:,1) = x0;

% Initial time and state
tk = tspan(1);
xk = x0;

for k = 1:N

    % Evaluate drift and diffusion
    f     = feval( driftFun    , tk, xk, u, d, p );
    sigma = feval( diffusionFun, tk, xk, u, d, p );
    
    % Advance with Euler-Maruyma
    tkp1 = tk + dt;
    xkp1 = xk + dt*f + sigma*dw(:, k);

    % Save solution
    T(k+1)    = tkp1;
    X(:, k+1) = xkp1;

    % Update point
    tk = tkp1;
    xk = xkp1;
    
end




