function [H_s_l, g_s_l, H_s_u, g_s_u] = getSoftOutputConstraints(N, model)
% Soft Output Constraints - Lower Bound 
W_hat_s_1_l =  1 * eye(model.su);
W_s_1_l = kron(eye(N), W_hat_s_1_l);

W_hat_s_2_l = 1 * eye(model.su);
W_s_2_l = kron(eye(N), W_hat_s_2_l);

H_s_l = W_s_2_l' * W_s_2_l;

e = ones(1, N * model.su);

g_s_l = W_s_1_l * e';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soft Output Constraints - Upper Bound 
W_hat_s_1_u =  1 * eye(model.su);
W_s_1_u = kron(eye(N), W_hat_s_1_u);

W_hat_s_2_u = 1 * eye(model.su);
W_s_2_u = kron(eye(N), W_hat_s_2_u);

H_s_u = W_s_2_u' * W_s_2_u;

g_s_u = W_s_1_u * e';


end

