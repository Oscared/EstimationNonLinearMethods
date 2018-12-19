%Iterative update [mean, Sigma] = update(mean_bar, Sigma_bar,
%z_bar, Q)

function [mean, Sigma] = update_i(mean_bar, Sigma_bar, z_bar, Q, dt)
    
    e = 1e-6;
    max_iter = 100;
    iter = 0;
    
    eta = mean_bar;
    eta_l = zeros(size(eta));
    
    while (norm(eta_l - eta,Inf)/norm(eta_l,Inf) >= e) && (iter < max_iter)
        
        if iter ~= 0
            eta = eta_l;
        end
        [h, H] = meas_model(eta, dt, Q);
        
        %Type 1
        %K = Sigma_bar*H'/(H*Sigma_bar*H' + Q);
        %eta_l = mean_bar + K*(z_bar - h' - H*(mean_bar - eta));
        %Type John
        Sigma_l = Sigma_bar - Sigma_bar*H'/(Q + H*Sigma_bar*H')*H*Sigma_bar;
        eta_l = eta + Sigma_l*(H'/Q*(z_bar - h') - Sigma_bar\(eta - mean_bar));
        %a = mod(mean(1,1)+pi,2*pi)-pi;
        %mean(1,1) = a;
        
        iter = iter + 1;
        
        %Sigma = (Sigma + Sigma')/2;
    end
    %disp(iter)
    mean = eta_l;
    %Type 1
    %Sigma = Sigma_bar - K*H*Sigma_bar;
    %Type John
    Sigma = Sigma_l;
end