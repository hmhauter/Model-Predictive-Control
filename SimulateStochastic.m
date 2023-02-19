function [y, z, x, d, t] = SimulateStochastic(F, scale, name)
    t0 = 0.0;
    tf = 20*60;
    Ts = 10;
    t = [t0:Ts:tf];     % time vector
    N = length(t);
    
    F1 = F(1);
    F2 = F(2);
    u = [repmat(F1,1,N); repmat(F2,1,N)];   % control 
    
    d = [randn(1,N); randn(1,N)];   % disturbances F3, F4
    parameters = GetParameters();
    
    P = eye(1);
    Pq = chol(P,'lower');

    % Process Noise 
    Q = scale * [10^2 0;0 20^2];
    Lq = chol(Q,'lower');
    w = Lq*randn(2,N);

    % Measurement Noise
    R = scale * eye(4);
    Lr = chol(R,'lower');
    v = Lr*randn(4,N);
    
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
        T = [T; Tk];
        X = [X; Xk];
    end
    
    k = N;
    y(:,k) = FourTankSystemSensor(x(:,k),parameters)+v(:,k); % Sensor function
    z(:,k) = FourTankSystemOutput(x(:,k),parameters);        % Output function
 
%     tl = strcat('Water Level in 4-Tank-System - ', name);
%     figure('Name', tl)
%     title(tl)
%     subplot(2,2,1)
%     plot(t,y(1, :))
%     title('Tank 1','interpreter','latex')
%     xlabel('Time [sec]','interpreter','latex')
%     ylabel('Water Level [cm]','interpreter','latex')
%     
%     subplot(2,2,2)
%     plot(t,y(2, :))
%     title('Tank 2','interpreter','latex')
%     xlabel('Time [sec]','interpreter','latex')
%     ylabel('Water Level [cm]','interpreter','latex')
%     
%     subplot(2,2,3)
%     plot(t,y(3, :))
%     title('Tank 3','interpreter','latex')
%     xlabel('Time [sec]','interpreter','latex')
%     ylabel('Water Level [cm]','interpreter','latex')
%     
%     subplot(2,2,4)
%     plot(t,y(4, :))
%     title('Tank 4','interpreter','latex')
%     xlabel('Time [sec]','interpreter','latex')
%     ylabel('Water Level [cm]','interpreter','latex')

end

