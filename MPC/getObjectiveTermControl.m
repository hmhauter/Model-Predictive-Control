function [H_u, g_u] = getObjectiveTermControl(N, model, U)
    % get objective functions
    % objective function for U
    % weight matrix

    % U 
    lambda_m = zeros(N*model.su^2, N*model.su^2);
    W_u = ones(model.su);
    W_hat_u = kron(eye(N), W_u);
    H_u = W_hat_u' * W_hat_u;
    g_u = - W_hat_u' * W_hat_u * repmat(U, N, 1);

end
