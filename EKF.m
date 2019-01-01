function [mu, Sigma] = EKF(mu, Sigma, dt, z, R, Q, constants)
% Function for doing the EKF one time, to play.m
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
%   mu               (2 or 4)x1     : new state efter the EKF process
%   Sigma     (2 or 4)x(2 or 4)     : new covariance matrix after the EKF

    % prediction
    g = process_model(mu, dt, constants);           % OBS not gravity
    G = process_model_jacobian(mu, dt, constants);
    
    mu_bar    = g;
    Sigma_bar = G*Sigma*G' + R;
    
    % measurement update
    h = measurement_model(mu_bar, constants);
    H = measurement_model_jacobian(mu_bar, constants);
    
    K = Sigma_bar*H'/(H*Sigma_bar*H' + Q);
    nu = z - h;
    
    sim_type = constants(11);
    meas_type = constants(12);  
        
    mu = mu_bar + K*nu;
    Sigma = Sigma_bar - K*H*Sigma_bar;

end
