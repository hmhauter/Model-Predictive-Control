function [x_k_k, w_k_k] = stationaryKalmanFilter(filter,x_prio_prio, w_prio_prio, u_prio_prio, y_k, d_k, model)
    
    x_k_prio = model.Ad * x_prio_prio + model.Bd * u_prio_prio + model.Ed *(w_prio_prio + d_k);
    y_k_prio = model.Cd * x_k_prio;
    e_k = y_k - y_k_prio;
    x_k_k = x_k_prio + filter.K_fx * e_k;
    w_k_k = filter.K_fw * e_k;
    
end
