function plot_errors(x_sim, x_ekf, x_iekf, x_ukf, Sigma_ekf, Sigma_iekf, Sigma_ukf, simulation_type)
% Function for plotting the errors from the filters in play.m
%
% INPUTS:     
%
%   x_sim            (2 or 4)xN     : simulated states
%   x_ekf            (2 or 4)xN     : estimated states with EKF
%   x_iekf           (2 or 4)xN     : estimated states with IEKF
%   x_ukf            (2 or 4)xN     : estimated states with UKF
%   simulation_type         int     : decides "2 or 4", 
%                                     0 = only the pendulum (then 2)
%                                     1 = cart and pendulum (then 4)
% OUTPUTS:
%
%   plots the errors of the filters compared the the simulation

% ----- Only pendulum (state is 2x1) -------------------------------------
if simulation_type == 0
   
    % ----- Error in angle of pendulum
    figure('Name', 'Error/Covariance plot')
    
    % EKF
    subplot(2,3,1)
    plot(mod((x_ekf(1,:) - x_sim(1,:)) + pi, 2*pi) - pi)
    %plot(x_ekf(1,:) - x_sim(1,:))
    title('Error EKF theta', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    %fprintf('MSE angle EKF: %f\n', mean((mod((x_ekf(1,:) - x_sim(1,:)) + pi, 2*pi) - pi)).^2)
    fprintf('MSE angle EKF: %f\n', mean((x_ekf(1,:) - x_sim(1,:)).^2))
    
    % IEKF
    subplot(2,3,2)
    plot(mod((x_iekf(1,:) - x_sim(1,:)) + pi, 2*pi) - pi)
    %plot(x_iekf(1,:) - x_sim(1,:))
    title('Error iEKF theta', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    %fprintf('MSE angle IEKF: %f\n', mean((mod((x_iekf(1,:) - x_sim(1,:)) + pi, 2*pi) - pi)).^2)
    fprintf('MSE angle IEKF: %f\n', mean((x_iekf(1,:) - x_sim(1,:)).^2))
    
    % UKF
    subplot(2,3,3)
    plot(mod((x_ukf(1,:) - x_sim(1,:)) + pi, 2*pi) - pi)
    %plot(x_ukf(1,:) - x_sim(1,:))
    title('Error UKF theta', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    %fprintf('MSE angle UKF: %f\n', mean((mod((x_ukf(1,:) - x_sim(1,:)) + pi, 2*pi) - pi)).^2)
    fprintf('MSE angle UKF: %f\n', mean((x_ukf(1,:) - x_sim(1,:)).^2))
    
    % ----- Covariance of angle 
    
    % EKF
    %figure('Name','Covariance angle')
    subplot(2,3,4)
    plot(reshape(Sigma_ekf(1,1,:), size(Sigma_ekf, 3), 1, 1))
    title('Variance angle EKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % IEKF
    subplot(2,3,5)
    plot(reshape(Sigma_iekf(1,1,:), size(Sigma_iekf, 3), 1, 1))
    title('Variance angle IEKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % UKF
    subplot(2,3,6)
    plot(reshape(Sigma_ukf(1,1,:), size(Sigma_ukf, 3), 1, 1))
    title('Variance angle UKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')

% ----- Cart and pendulum (state is 4x1) ---------------------------------    
elseif simulation_type == 1
    
    % ----- Error in cart position
    figure('Name', 'Error plot')
    
    % EKF
    subplot(2,3,1)
    plot(((x_ekf(1,:) - x_sim(1,:))))
    title({'Error', 'x EKF'}, 'interpreter', 'latex')
    %title('Error x EKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    fprintf('MSE pos cart EKF: %f\n', mean((x_ekf(1,:) - x_sim(1,:)).^2))
    
    % IEKF
    subplot(2,3,2)
    plot((x_iekf(1,:) - x_sim(1,:)))
    title({'Error', 'x IEKF'}, 'interpreter', 'latex')
    %title('Error x IEKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    fprintf('MSE pos cart IEKF: %f\n', mean((x_iekf(1,:) - x_sim(1,:)).^2))
    
    % UKF
    subplot(2,3,3)
    plot((x_ukf(1,:) - x_sim(1,:)))
    title({'Error', 'x UKF'}, 'interpreter', 'latex')
    %title('Error x UKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    fprintf('MSE pos cart UKF: %f\n', mean((x_ukf(1,:) - x_sim(1,:)).^2))

    % ----- Error in angle of pendulum
    
    % EKF
    subplot(2,3,4)
    plot((mod((x_ekf(3,:) - x_sim(3,:)) + pi, 2*pi) - pi))
    %plot(x_ekf(3,:) - x_sim(3,:))
    title('Error $\theta$ EKF' ,'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    %fprintf('MSE angle EKF: %f\n', mean((mod((x_ekf(3,:) - x_sim(3,:)) + pi, 2*pi) - pi)).^2)
    fprintf('MSE angle EKF: %f\n', mean((x_ekf(3,:) - x_sim(3,:)).^2))
    
    % IEKF
    subplot(2,3,5)
    plot((mod((x_iekf(3,:) - x_sim(3,:)) + pi, 2*pi) - pi))
    %plot(x_iekf(3,:) - x_sim(3,:))
    title('Error $\theta$ IEKF' ,'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    %fprintf('MSE angle IEKF: %f\n', mean((mod((x_iekf(3,:) - x_sim(3,:)) + pi, 2*pi) - pi)).^2)
    fprintf('MSE angle IEKF: %f\n', mean((x_iekf(3,:) - x_sim(3,:)).^2))
    
    %UKF
    subplot(2,3,6)
    plot((mod((x_ukf(3,:) - x_sim(3,:)) + pi, 2*pi) - pi))
    %plot(x_ukf(3,:) - x_sim(3,:))
    title('Error $\theta$ UKF' ,'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    %fprintf('MSE angle UKF: %f\n', mean((mod((x_ukf(3,:) - x_sim(3,:)) + pi, 2*pi) - pi)).^2)
    fprintf('MSE angle UKF: %f\n', mean((x_ukf(3,:) - x_sim(3,:)).^2))
    
    % ----- Covariance of cart position 
    figure('Name', 'Covariance of cart position')
    
    % EKF
    subplot(2,3,1)
    plot(reshape(Sigma_ekf(1,1,:), size(Sigma_ekf, 3), 1, 1))
    title({'Variance', 'x EKF'}, 'interpreter', 'latex')
    %title('Variance x EKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % IEKF
    subplot(2,3,2)
    plot(reshape(Sigma_iekf(1,1,:), size(Sigma_iekf, 3), 1, 1))
    title({'Variance', 'x IEKF'}, 'interpreter', 'latex')
    %title('Variance x IEKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % UKF
    subplot(2,3,3)
    plot(reshape(Sigma_ukf(1,1,:), size(Sigma_ukf, 3), 1, 1))
    title({'Variance', 'x UKF'}, 'interpreter', 'latex')
    %title('Variance x UKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % ----- Covariance of angle of the pendulum
    
    % EKF
    subplot(2,3,4)
    plot(reshape(Sigma_ekf(3,3,:), size(Sigma_ekf, 3), 1, 1))
    title({'Variance', '$\theta$ EKF'}, 'interpreter', 'latex')
    %title('Variance $\theta$ EKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % IEKF
    subplot(2,3,5)
    plot(reshape(Sigma_iekf(3,3,:), size(Sigma_iekf, 3), 1, 1))
    title({'Variance', '$\theta$ IEKF'}, 'interpreter', 'latex')
    %title('Variance $\theta$ IEKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
    % UKF
    subplot(2,3,6)
    plot(reshape(Sigma_ukf(3,3,:), size(Sigma_ukf, 3), 1, 1))
    title({'Variance', '$\theta$ UKF'}, 'interpreter', 'latex')
    %title('Variance $\theta$ UKF', 'interpreter', 'latex')
    xlabel('Step [nr]','interpreter','latex')
    
end