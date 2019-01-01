% _____ Inverted pendulum simulation and estimation ______________________

close all
clear

% ----- Type of simulation -----------------------------------------------
%sim_type = 0;  % pendulum, 2D state
sim_type = 1;   % pendulum + cart, 4D state

% ----- Type of measurement ----------------------------------------------
%meas_type = 0;  % angle, pos (linear)
meas_type = 1;  % x, y for angle, pos (non-linear)

% ----- What plots? ------------------------------------------------------
plot_type = 0;   % simulation 
%plot_type = 1;   % simulation + errors
%plot_type = 2;   % errors


% ----- Constants --------------------------------------------------------
% time/step
Tmax = 100;
dt = 0.05;

% constants = [ g, l, friction,  F,  m,   M,  w,   h,      r,   lim, sim_type  meas_type]
%               1  2      3      4   5    6   7    8       9     10     11
constants = [9.82, 1,   0.01,   0,  1,   5, 0.5, 0.5/3, 0.5/3, 1.5, sim_type, meas_type];

% simulation
dt_sim = dt;
N_sim = Tmax/dt_sim;

% estimation
dt_est = dt*4;
N_est = Tmax/dt_est;

% relation of time steps
rel = dt_est/dt_sim;

% ----- Variances --------------------------------------------------------
R = 0;
Q = 0;

if sim_type == 0                           % linear measurement
    
    R = diag([0.01, 0.01]);                       % process
    
    if meas_type == 0                          % only pendulum
        
        Q = 0.1;                                  % measurement
        
    elseif meas_type == 1                      % with cart
        
       Q = diag([0.1, 0.1]);                     % measurement
       
    end    
    
elseif sim_type == 1                       % non-linear measurement
    
    R = diag([0.1, 0.1, 0.1, 0.1]);            % process
    
    if meas_type == 0                          % only pendulum
        
        Q = diag([0.01, 0.01]);                      % measurement
        
    elseif meas_type == 1                      % with cart
        
        Q = diag([0.01, 0.01, 0.01]);           % measurement
        
    end
    
end    

% ----- Initialization ---------------------------------------------------
x0 = 0;
Sigma0 = 0;
sz_x = 0;

% only pendulum
if sim_type == 0
    
    x0 = [1; 0];
    Sigma0 = diag([0.01, 0.01]);
    sz_x = size(x0,1);

% pendulum with cart    
elseif sim_type == 1
   
    x0 = [0; 0; 1; 0];
    Sigma0 = diag([0.01, 0.001, 0.01, 0.001]);
    sz_x = size(x0,1);
    
end

% ----- Simulation -------------------------------------------------------
x_sim = zeros(sz_x, N_sim);
x_sim(:,1) = x0;

for i = 2:N_sim
    
    x_sim(:,i) = process_model(x_sim(:,i-1), dt_sim, constants);
          
end

% ----- Measurements -----------------------------------------------------
x_sim_short = x_sim(:,1:rel:end);
x_meas = 0;
Q_noise = 1*Q;

% linear (measure angle + position of cart)
if meas_type == 0
   
    x_meas = (x_sim_short(1:2:end,:) + randn(size(Q, 1), N_est).*sqrt(diag(Q_noise)));

% non-linear (measure x,y from angle + position of cart)    
elseif meas_type == 1
    
    sim = x_sim_short;
    l = constants(2);
    
    % only pendulum
    if sim_type == 0
    
        x_meas = [l*cos(sim(1,:)) + randn(1, N_est)*sqrt(Q_noise(1,1));
                  l*sin(sim(1,:)) + randn(1, N_est)*sqrt(Q_noise(1,1))];
           
    % with cart
    elseif sim_type == 1
        
        x_meas = [sim(1,:) + randn(1, N_est)*sqrt(Q(1,1));
                  sim(1,:) + l*cos(sim(3,:)) + randn(1, N_est)*sqrt(Q_noise(2,2));
                  sim(1,:) + l*sin(sim(3,:)) + randn(1, N_est)*sqrt(Q_noise(2,2))];
    
    end
end    

% ----- Estimation -------------------------------------------------------

% EKF
x_ekf            = zeros(sz_x, N_est);
x_ekf(:,1)       = x0;
Sigma_ekf        = zeros(sz_x,sz_x,N_est);
Sigma_ekf(:,:,1) = Sigma0;

% IEKF
x_iekf            = zeros(sz_x, N_est);
x_iekf(:,1)       = x0;
Sigma_iekf        = zeros(sz_x,sz_x,N_est);
Sigma_iekf(:,:,1) = Sigma0;

% UKF
x_ukf            = zeros(sz_x, N_est);
x_ukf(:,1)       = x0;
Sigma_ukf        = zeros(sz_x,sz_x,N_est);
Sigma_ukf(:,:,1) = Sigma0;

for i = 2:N_est
    
    % measurement
    z = x_meas(:,i);
    
    % EKF
    [x_ekf(:,i), Sigma_ekf(:,:,i)] = EKF(x_ekf(:,i-1), Sigma_ekf(:,:,i-1), ...
        dt_est, z, R, Q, constants);
  
    % IEKF
    [x_iekf(:,i), Sigma_iekf(:,:,i)] = IEKF(x_iekf(:,i-1), Sigma_iekf(:,:,i-1), ...
        dt_est, z, R, Q, constants);
    
    % UKF
    [x_ukf(:,i), Sigma_ukf(:,:,i)] = UKF(x_ukf(:,i-1), Sigma_ukf(:,:,i-1), ...
        dt_est, z, R, Q, constants);
    
end

% ----- Plot the simulation with estimates -------------------------------
if plot_type < 2

    figure(1)

    for i = 1:N_est   

       plot_system(x_sim_short(:,i), i, constants);    % simulation
       plot_estimation(x_ekf(:,i),'r', constants);     % EKF
       plot_estimation(x_iekf(:,i),'m', constants);    % IEKF
       plot_estimation(x_ukf(:,i),'g', constants);     % UKF
       
       title(sprintf('Step %d of %d', i, N_est))
       L(1) = plot(NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'linewidth', 4);
       L(2) = plot(NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'linewidth', 4);
       L(3) = plot(NaN, NaN, 'go', 'MarkerFaceColor', 'g', 'linewidth', 4);
       legend(L, {'EKF', 'IEKF', 'UKF'})
       pause(dt_est)
       %pause(0.0001)

    end
    
end   

% ----- Plot errors ------------------------------------------------------
if plot_type > 0
    
    plot_errors(x_sim_short, x_ekf, x_iekf, x_ukf, Sigma_ekf, Sigma_iekf, Sigma_ukf, constants(11));

end

