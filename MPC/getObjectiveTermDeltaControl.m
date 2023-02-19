function [H_du, g_du, lambda_m] = getObjectiveTermDeltaControl(N, model, u_prev)
    % get objective functions
    % objective function for z
    % setpoints r
    % weight matrix

    % rate of movement delta U 
    lambda_m = zeros(N*model.su^2, N*model.su^2);
%     W_du = 0.1 * ones(model.su);
    W_du =  eye(model.su);
    W_hat_du = kron(eye(N), W_du);

    helper_1 = ones(1, model.su*N);
    helper_2 = -1*ones(1, model.su*(N-1));
    lambda_m = diag(helper_1) + diag(helper_2, -model.su);
    eye_0 = repmat(eye(model.su),N, 1);

    H_du = (W_hat_du*lambda_m)'*(W_hat_du*lambda_m);

    g_du  = -(W_hat_du* lambda_m)'*W_hat_du*eye_0*u_prev; %repmat(u_prev, N ,1);
    
    % M_z = -(W_hat_z*kalmanFilter.gamma)'*W_hat_z;
end