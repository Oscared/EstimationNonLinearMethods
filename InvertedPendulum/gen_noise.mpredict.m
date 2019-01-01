%Predict [mean_bar, Sigma_bar] = predict(x, Sigma, dt, friction, R, g, l)
%Gives prediction of movement

function [mean_bar, Sigma_bar] = predict(x, Sigma, dt, friction, R, g, l)
    
    mean_bar = x + [dt*x(2,1); (g/l*sin(x(1,1))*dt - friction*x(2,1))];% + sqrt(diag(R)).*ones(2,1);
    %a = mod(mean_bar(1,1)+pi,2*pi)-pi;
    %mean_bar(1,1) = a;
    G = [1, dt; g*dt/l*cos(x(1,1)), 1-friction];
    Sigma_bar = G*Sigma*G' + R;
        
end
