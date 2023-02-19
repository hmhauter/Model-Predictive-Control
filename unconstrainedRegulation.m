function u_opt = unconstrainedRegulation(model, u_prev, N, ...
    omega_x, omega_w, gamma, z_est, x_est, y_est, x_k, w_k)
% N: number of predicion steps 

    % [z_est, x_est, y_est]  = kalmanPrediction(dx, dy, du, w, v, N, model, kalman);

    [H_z, g_z] = getObjectiveTermOutput(omega_x, omega_w, gamma, N, model, x_k, w_k);

    [H_du, g_du] = getObjectiveTermDeltaControl(N, model, u_prev);

    H = H_z + H_du;
    g = g_z + g_du;     % DIMENSION MISTAKE HERE 

    [u_opt,fval,exitflag,output,lambda] = quadprog(H,g)
    disp(u_opt)
    disp(fval)
    disp(output)

%     u_opt = regulator(H,g,l,u,A,bl,bu,xinit)

end

function [dx, dy] = simulateOneStep(model, u_opt, x, w, F, v)
    dx = model.Ad*x+model.Bd*u_opt + model.Ed*(w+F(3:4)');
    dy = model.Cd*x + v;

end

function [H_z, g_z] = getObjectiveTermOutput(omega_x, omega_w, gamma, N, model, x_pred, w_pred)
    % get objective functions
    % objective function for z
    % setpoints r
    % weight matrix 
    W_z = ones(N);
    W_hat_z = kron(eye(N), W_z);
    R_k =  [];
    % model.Rk;  % setpoints
    % repeat R_k N times 
    for i = 1:N
       R_k = [R_k; model.r']; 
    end
    
    b_k = omega_x * x_pred + omega_w * w_pred(1:2, :);
    
    c_k = R_k - b_k;

    H_z = (W_hat_z * gamma)'*(W_hat_z*gamma);
    g_z  = -(W_hat_z*gamma)'*W_hat_z*c_k;
    
    % M_z = -(W_hat_z*kalmanFilter.gamma)'*W_hat_z;
end


function [H_du, g_du] = getObjectiveTermDeltaControl(N, model, u_prev)
    % get objective functions
    % objective function for z
    % setpoints r
    % weight matrix

    % rate of movement delta U 
    lambda_m = zeros(N*model.su^2, N*model.su^2);
    W_du = ones(model.su);
    W_hat_du = kron(eye(N), W_du);

    helper_1 = ones(1, model.su*N);
    helper_2 = -1*ones(1, model.su*(N-1));
    lambda_m = diag(helper_1) + diag(helper_2, -model.su);
    eye_0 = repmat(eye(model.su),N, 1);

    H_du = (W_hat_du*lambda_m)'*(W_hat_du*lambda_m);

    g_du  = -(W_hat_du* lambda_m)'*W_hat_du*eye_0*u_prev; %repmat(u_prev, N ,1);
    
    % M_z = -(W_hat_z*kalmanFilter.gamma)'*W_hat_z;
end

function [z_est, x_est, y_est] = kalmanPrediction(dx, dy, du, w, v, N, model, kalman)

    for k = 1:N
        dx(:,k+1) = model.Ad*dx(:,k)+model.Bd*du + model.Ed*(w(1:2,k)+F(3:4)');
        dy(:,k+1) = model.Cd*dx(:,k+1) + v(:,k+1);
    
        e = dy(:,k+1) - model.Cd*dx(:,k+1);
        
        x_est(:, k) = dx(:,k+1) + kalman.K_fx*e;
        y_est(:,k) = model.Cd*x_est;
        z_est(:, k) = model.Czd * x_est;
    end

end

% function [omega_x, omega_w, gamma] = setupMatrices(model, N)
%     % N is the number of prediction steps
%     % allocate memory
%     omega_x = zeros(N*model.sz, model.sz);
%     omega_w = zeros(N*model.sz, model.su);
%     omega_h = zeros(N*model.sz, model.su);
% 
%     H = zeros(N*model.sz, N*model.su);
% 
%     omega_x(1:model.sz, :) = model.Cd*model.Ad;
%     omega_w(1:model.sz, :) = model.Cd*model.Ed;
%     omega_h(1:model.sz, :) = model.Cd*model.Bd;
% 
%      
%     % rows
%     for j = 2:N
%         omega_x((j-1)*model.sz+1:j*model.sz, :) = omega_x((j-2)*model.sz+1:(j-1)*model.sz, 1:model.sz)*model.Ad;
%         omega_w((j-1)*model.sz+1:j*model.sz, :) = omega_x((j-2)*model.sz+1:(j-1)*model.sz, 1:model.sz)*model.Ed;
%         omega_h((j-1)*model.sz+1:j*model.sz, :) = omega_x((j-2)*model.sz+1:(j-1)*model.sz, 1:model.sz)*model.Bd;
%     end
% 
%     % post processing for H matrix 
%     % columns
%     H(1:model.sz*model.sz, 1:model.su) = omega_h;
%     for j = 2:N
%         H((j-1)*model.sz+1:N*model.sz, (j-1)*model.su+1:j*model.su) = omega_h(1:(N-j+1)*model.sz,:);
%     end
% 
%     gamma = H;
% 
% end