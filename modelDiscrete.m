function [Model] = modelDiscrete(Ad, Bd, Cd, Czd, Ed, yss)
    Model.Ad = Ad;
    Model.Bd = Bd;
    Model.Cd = Cd;
    Model.Czd = Czd;
    Model.Ed = Ed;

    % Process Noise 
    Q = [10^2 0;0 20^2];
    Lq = chol(Q,'lower');

    Model.Qk = Q;
    
    % Measurement Noise
    R = 25 * eye(4);
    Lr = chol(R,'lower');

    Model.Rk = R;


    
    Model.Sk = zeros(2,4);

    Model.su = 2;
    Model.sz = 2;
    Model.sx = 4;

    Model.r = [100, 110, 40, 40]-yss';


end


