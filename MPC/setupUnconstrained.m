function [Lx, Lr, Lu, Ld] = setupUnconstrained(model,S, Q, N)
    % N is the number of prediction steps
    % allocate memory
    omega_x = zeros(N*model.su, model.sx);
    omega_w = zeros(N*model.su, model.su);
    omega_h = zeros(N*model.su, model.su);

    H = zeros(N*model.su, N*model.su);

    omega_x(1:model.su, :) = model.Czd*model.Ad;
    omega_w(1:model.su, :) = model.Czd*model.Ed;
    omega_h(1:model.su, :) = model.Czd*model.Bd;

     
    % rows
    for j = 2:N
        omega_x((j-1)*model.su+1:j*model.su, :) = omega_x((j-2)*model.su+1:(j-1)*model.su, 1:model.sx)*model.Ad;
        omega_w((j-1)*model.su+1:j*model.su, :) = omega_x((j-2)*model.su+1:(j-1)*model.su, 1:model.sx)*model.Ed;
        omega_h((j-1)*model.su+1:j*model.su, :) = omega_x((j-2)*model.su+1:(j-1)*model.su, 1:model.sx)*model.Bd;
    end

    % post processing for H matrix 
    % columns
    % Cd * Bd
    H(:, 1:model.su) = omega_h(:, :);
    Hd(:, 1:model.su) = omega_w(:, :);
    for j = 2:N
        H((j-1)*model.su+1:N*model.su, (j-1)*model.su+1:j*model.su) = omega_h(1:(N-j+1)*model.su,:);
        Hd((j-1)*model.su+1:N*model.su, (j-1)*model.su+1:j*model.su) = omega_w(1:(N-j+1)*model.su,:);
    end

    gamma = H;
    gamma_w = Hd;


    lambda_m = 2*diag(ones(1,N*model.su), 0) - diag(ones(1,(N-1)*model.su), -model.su) - diag(ones(1,(N-1)*model.su), model.su);
    
    % get the gains for the MPC
    I0 = eye(N*model.su,model.su);
    I0 = I0';


    Q_z = kron(eye(N,N),Q);
%     S_u = kron(eye(model.su,model.su),S);
    S_u = kron(eye(N/model.su,N/model.su),S);
    Hp = gamma' * Q_z * gamma + lambda_m;
    M_x0 = gamma' * Q_z * omega_x;
    M_r = gamma' * Q_z;
    M_d = gamma'*Q_z*gamma_w;
    M_u = S_u;
    H_inv = inv(Hp);
    Lx = -H_inv * M_x0;
    Lr = H_inv * M_r;
    Lu = -H_inv * M_u;
    Ld = -H_inv * M_d;
    
    Lx = I0 * Lx;
    Lr = I0 * Lr;
    Lu = I0 * Lu;
    Lu = Lu(1:2, 1:2);
    Ld = I0 * Ld;
    Ld = Ld(1:2, 1:2);
end

