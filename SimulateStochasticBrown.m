function [y, z, t] = SimulateStochasticBrown(F)
    t0 = 0.0;
    tfin = 20*60;
    Ts = 10;
    t = [t0:Ts:tfin];     % time vector
    N = length(t);
    
    sz = 4;
    su = 2;
    
    % Realization of Ns Standard Brownian Motions
    Ns = 4;            % Number of realizations
    seed = 100;         % Seed for reproducibility
    rng(seed);
    dV = sqrt(Ts)*randn(Ns,N);
    v = cumsum(dV,2);
    
    F1 = F(1);
    F2 = F(2);
    u = [repmat(F1,1,N); repmat(F2,1,N)];   % control 
    
    d = [randn(1,N); randn(1,N)];   % disturbances F3, F4
    r = [13, 10, 4, 4];
    parameters = GetParameters();
    
    P = eye(1);
    Pq = chol(P,'lower');
    x0 = Pq*randn(4, 1);
    
    % Process Noise 
    Q = [1^2 0;0 1^2];
    Lq = chol(Q,'lower');
    w = Lq*randn(2,N);

    nx = 4; nu = 2; ny = 4; nz = 2;
    x = zeros(nx,N);
    y = zeros(ny,N);
    
    
    z = zeros(nz,N);
    X = zeros(0,nx);
    T = zeros(0,1);
    step_size = 0.01;
    x(:,1) = [0;0;0;0];
    
    for k = 1:N-1
        y(:,k) = FourTankSystemSensor(x(:,k),parameters)+v(:,k); % Sensor function WITH noise 
        z(:,k) = FourTankSystemOutput(x(:,k),parameters);
        % Output function
        disturbance_with_noise = d(:,k)+w(:,k);
        [Tk,Xk] = ode15s(@ModifiedFourTankSystem,[t(k):step_size:t(k+1)],x(:,k),[], ...
            u(:,k),disturbance_with_noise,parameters);
        x(:,k+1) = Xk(end,:)';
        % T = [T; Tk];
        X = [X; Xk];
    end
    
    k = N;
    y(:,k) = FourTankSystemSensor(x(:,k),parameters)+v(:,k); % Sensor function
    z(:,k) = FourTankSystemOutput(x(:,k),parameters);        % Output function


end

