%Measurement model [h, H] = meas_model(mean_bar, dt, Q)

function [h, H] = meas_model(mean_bar, dt, Q)
    l = 1;
    %h = mean_bar(1,1) + mean_bar(2,1)*dt; % + sqrt(Q)*randn;
    %H = [1, dt];
    %Linear model with state of theta and theta dot
    %h = mean_bar(1,1);
    %H = [1,0];
    %Non-linear measurement model
    h = [l*sin(mean_bar(1,1)), l*cos(mean_bar(1,1))];
    H = [l*cos(mean_bar(1,1)), 0; -l*sin(mean_bar(1,1)), 0];
end