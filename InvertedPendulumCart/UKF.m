function [mu, Sigma] = UKF(mu, Sigma, dt, z, R, Q, constants)
% Function for doing the UKF one time, to play.m
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
%   mu               (2 or 4)x1     : new state efter the UKF process
%   Sigma     (2 or 4)x(2 or 4)     : new covariance matrix after the UKF
    
    % ----- Number of dimensions -----------------------------------------
    n = size(mu,1);
    
    % ----- Dimension of measurement
    m = size(z, 1);

    % ----- UKF constants ------------------------------------------------
    alpha  = 1e-3; 
    beta   = 0;
    kappa  = 2;
    
    lambda = alpha^2*(n + kappa) - n;
    gamma  = sqrt(n + lambda);
    
    % ----- Weights ------------------------------------------------------
    wm0 = lambda/(n + lambda);
    wc0 = wm0 + (1 - alpha^2 + beta);
    
    wmci = 1/(2*(n + lambda));
    
    Wm = [wm0, wmci*ones(1,2*n)];
    Wc = [wc0, wmci*ones(1,2*n)];
    
    % ----- Sigma points -------------------------------------------------
    X = sigma_points(mu, Sigma, gamma);
    
    % ----- Unscented transform of process --------------------------------
    [mu_bar, Sigma_bar, ~] = unscented_transform(X, Wm, Wc, n, n, R, dt, constants, 0);
    
    % ----- Sigma points of prediction -----------------------------------
    mu_bars = sigma_points(mu_bar, Sigma_bar, gamma);
    mu_diff = mu_bars - mu_bar(:,ones(1,2*n+1));
    
    % ----- Unscented transform of measurements ---------------------------
    [z_hat, S, z_diff] = unscented_transform(mu_bars, Wm, Wc, m, n, Q, dt, constants, 1);
    
    % ----- Sigma_xz (crosscovariance) -----------------------------------
    Sigma_xz = mu_diff*diag(Wc)*z_diff';

    
    % ----- K and final results ------------------------------------------
    K = Sigma_xz/S;
    nu = z - z_hat;
    
    sim_type = constants(11);
    meas_type = constants(12);
    % bound between -pi and pi
%     if meas_type == 0
%         if sim_type == 0
%             nu = mod(nu + pi, 2*pi) - pi;
%         elseif sim_type == 1
%             nu(2) = mod(nu(2) + pi, 2*pi) - pi;
%         end
%     end
    
    mu = mu_bar + K*nu;
    Sigma = Sigma_bar - K*Sigma_xz';

end

function [x, S, x_diff] = unscented_transform(pts, Wm, Wc, n, N, P, dt, constants, stage)
% Function for doing the unscented transform to the UKF function
%
% INPUTS:     
%
%   pts              (2 or 4)xn     : sigma points
%   Wm                 1x(2n+1)     : weights for the mean
%
% OUTPUTS:
%
%   mu               (2 or 4)x1     : new state efter the UKF process
%   Sigma     (2 or 4)x(2 or 4)     : new covariance matrix after the UKF
    
    % type: 1 or 0, 1 = process model, 0 = measurement model

    % X: sigma points that went though the process / measurement model 
    % x: z_hat, mu_bar   -- weighted sum of X
    % S: Sigma_bar/S : covariance of prediction/measurement
    % X_diff: difference between X and x

    x = zeros(n,1);
    X = zeros(n, 2*N+1);
    
    for k = 1:2*N+1
        
       if stage == 0
           
           X(:,k) = process_model(pts(:,k), dt, constants);
       
       elseif stage == 1
          
           X(:,k) = measurement_model(pts(:,k), constants);
       
       end
       
       % mean
       x = x + Wm(k)*X(:,k);
        
    end
    
     x_diff = X - x(:,ones(1,2*N+1));
    
    % variance
    S = x_diff*diag(Wc)*x_diff' + P;
    

end


function X = sigma_points(mu, Sigma, gamma)
% Function for doing generating the sigma points, to UKF
%
% INPUTS:     
%
%   mu               (2 or 4)x1     : old state
%   Sigma     (2 or 4)x(2 or 4)     : old covariance
%   gammma                float     : related to the spread of the sigma
%                                     points
%
% OUTPUTS:
%
%   X           (2 or 4)x(2n+1)     : sigma points

    A = gamma*chol(0.5*(Sigma + Sigma'))';
    Y = mu(:,ones(1, size(mu, 1)));
    X = [mu, Y+A, Y-A];

end    
