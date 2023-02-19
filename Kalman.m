function [Kalman] = Kalman(Ad, Cd)
% set up matrixes for Kalman Filter
    G = eye(4);
    R = eye(4); 
    Q = 5^2 * eye(4);
    S = eye(4);
    % R_ww = R(1,1); R_wv = R(1,2); R_vw = R(2,1); R_vv = R(2,2);
    % P = dare(Ad',Cd',G*Q*G',R_vv,G*R_wv);
    P = dare(Ad',Cd',G*Q*G',R,G*S);
    
    R_e = Cd * P * Cd' + R;
    K_fx = P * Cd' * pinv(R_e);
    
    K_fw = S * pinv(R_e);

    Kalman.K_fx = K_fx;
    Kalman.K_fw = K_fw;
    Kalman.S = S;
    Kalman.P = P;

end

