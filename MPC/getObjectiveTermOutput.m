function [H_z, g_z] = getObjectiveTermOutput(omega_x, omega_w, gamma, N, r_set, x_pred, w_pred)
    % get objective functions
    % objective function for z
    % setpoints r
    % weight matrix 
    W_z =  1000 * eye(N);% ones(N);
    W_hat_z = kron(eye(4), W_z);
    R_k =  [];
    % model.Rk;  % setpoints
    % repeat R_k N times 
    
    b_k = omega_x * x_pred + omega_w * w_pred(1:2, :);
    
    c_k = r_set' - b_k;

    H_z = (W_hat_z * gamma)'*(W_hat_z*gamma);
    g_z  = -(W_hat_z*gamma)'*W_hat_z*c_k;
    
    % M_z = -(W_hat_z*kalmanFilter.gamma)'*W_hat_z;
end