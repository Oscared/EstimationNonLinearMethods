function [mu, Sigma] = IEKF(mu, Sigma, dt, z, R, Q, constants)
% Function for doing the IEKF one time, to play.m
%
% INPUTS:     
%
%   mu               (2 or 4)x1     : old state
%   Sigma     (2 or 4)x(2 or 4)     : old covariance
%   dt                    float     : time step
%   z              (1,2 or 3)x1     : measurement
%   R         (2 or 4)x(2 or 4)     : process noise covariance
%   Q     (1,2 or 3)x(1,2 or 3)     : measurement noise covariance
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   mu               (2 or 4)x1     : new state efter the IEKF process
%   Sigma     (2 or 4)x(2 or 4)     : new covariance matrix after the IEKF

    % prediction
    g = process_model(mu, dt, constants);           % OBS not gravity
    G = process_model_jacobian(mu, dt, constants);
    
    mu_bar    = g;
    Sigma_bar = G*Sigma*G' + R;
    
    % ----- Iterative part -----------------------------------------------
    % maximum number of iterations
    max_iter = 100;

    % tolerance
    e = 1e-10;
    
    mu = mu_bar;
    mu_last = mu;
    
    iter = 0;
    while (iter == 0) || ((norm(mu - mu_last, Inf) > e) && (iter < max_iter))
       
        mu_last = mu;
        H = measurement_model_jacobian(mu, constants);
        Sigma = Sigma_bar - Sigma_bar*H'/(H*Sigma_bar*H' + Q)*H*Sigma_bar;
        nu = z - measurement_model(mu, constants);
        
        sim_type = constants(11);
        meas_type = constants(12);
        
        mu = mu + Sigma*H'/Q*(nu) - Sigma/Sigma_bar*(mu - mu_bar);
        iter = iter + 1;
        
    end
    
    H = measurement_model_jacobian(mu, constants);
    Sigma = Sigma_bar - Sigma_bar*H'/(H*Sigma_bar*H' + Q)*H*Sigma_bar;
    
end
