%Update step for the iEKF [mean, Sigma] = update(mean_bar, Sigma_bar,
%z_bar, Q)

function [mean, Sigma] = update_(mean_bar, Sigma_bar, z_bar, Q, h, H)
    
    K = Sigma_bar*H'/(H*Sigma_bar*H' + Q);
    mean = mean_bar + K*(z_bar - h);
    %a = mod(mean(1,1)+pi,2*pi)-pi;
    %mean(1,1) = a;
    Sigma = Sigma_bar - K*H*Sigma_bar;
    
    %Sigma = (Sigma + Sigma')/2;
    
end
