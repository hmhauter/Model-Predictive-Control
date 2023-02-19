function [z, x, R, x_hat_post_k, P_post_k] = dynamicKalmanFilter(x_hat_prio,P_prio, y_k, u_k, d_k, model, N) 
    
    % compute innovation 
    y_hat_prio = model.Cd*x_hat_prio;
    e_k = y_k - y_hat_prio;
    R_e_k = model.Cd * P_prio * model.Cd' + model.Rk;
    R_e_k_inv = inv(R_e_k);

    % compute filtered state and filtered process noise 
    K_fx_k = P_prio*model.Cd'*R_e_k_inv;
    K_fw_k = model.Sk*R_e_k_inv;

    x_hat_k = x_hat_prio + K_fx_k * e_k;
    w_hat_k = K_fw_k * e_k;
    P_k_k = P_prio - K_fx_k * R_e_k*K_fx_k';
    Q_k_k = model.Qk - K_fw_k * R_e_k * K_fw_k';
    
    % state predictions
    x_hat_post_k = model.Ad * x_hat_k + model.Bd * u_k + model.Ed * (d_k + w_hat_k);
    P_post_k = model.Ad * P_k_k * model.Ad' + model.Ed * model.Qk * model.Ed';
    
%     P = cell(N, 1);
%     R = cell(N, 1);
%     x = zeros(4,N);    
%     z = zeros(4,N);
% 
%     x(:, 1) = x_hat_post_k;
%     P{1} = P_post_k;
% 
% 
%     j = 1;
%     for i = 2:N-1
%         x(:, i) = model.Ad * x(:, i-1);
%         P{i} = model.Ad * P{i-1} * model.Ad' + model.Ed * model.Qk * model.Ed';
%         
%         z(:, j) = model.Cd * x(:, j);
%         R{j} = model.Cd * P{j} * model.Cd';
% 
%         j = j+1;
%     end
     
    
    

end
