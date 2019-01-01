%EKF test
close all
%Initialization 
N = 40;

%If plot moving set to 1
plot_setting = 0;

%If do uiKF set to 1
uikf_setting = 0;

% --- Process ---
R = diag([0.01, 0.001]);
dt = 0.1;
g = 9.82;
l = 1;
friction = 0.01;
N_ekf = N/dt;

% --- Measurement ---
Q = 0.1;

% --- Simulation ---
dt_sim = 0.01;
N_sim = N/dt_sim;

% --- Start values ---
x_start =[0; 0.1];
Sigma_bar = diag([0.01 ,0.001]);


%% Simulation
th = zeros(N_sim,2);
th(1,2) = 0.1;

for i=1:N_sim-1
    a = th(i,1) + th(i,2)*dt_sim;
    %a = mod(a+pi,2*pi)-pi;
    th(i+1,1) = a;
    th(i+1,2) = th(i,2) + g/l*sin(th(i,1))*dt_sim - friction*th(i,2);
end    

%% Simulation like Mikael

th_m = zeros(N_sim,2);
th_m(1,2) = x_start(2,1);

for i=2:N_sim
    th_m(i,:) =  th_m(i-1,:) + [dt_sim*th_m(i-1,2), (g/l*sin(th_m(i-1,1))*dt_sim - friction*th_m(i-1,2))];
end

%th = th_m;

%% Testing other non-linear model
%Measuring x and y instead of theta and theta dot
Q = 0.1.*eye(2);
th_linear = th;
th = th(1:dt/dt_sim:end,:);
%Generate noise function
noise = gen_noise(th,0,Q);
th = th + noise;
%Linear measurement 
%z_bar = th(:,1);
%Non-linear measurement
z_bar = [l.*sin(th(:,1)), l.*cos(th(:,1))];
%Noise after transform
%z_bar = [l.*sin(th(:,1)), l.*cos(th(:,1))] + noise;

%Model Q
Q = 0.1*eye(2);

%Giving noise to parameters to change prediciton
new = [friction, g, l] + 0.5.*[friction, g, l].*randn(1,3);
friction = new(1); g = new(2); l = new(3);
%% EKF
x = zeros(2, N_ekf);
x(:,1) = x_start;
%Linear measuremeent with theta and theta dot
%z_bar = th(:,1) + sqrt(Q).*randn(size(th(:,1)));
%Measurment non-linear as x and y
%z_bar = th(1:dt/dt_sim:end,:) + (sqrt(diag(Q))'.*randn(N_ekf,2));

P_e = zeros(2,2,N_ekf);
P_e(:,:,1) = Sigma_bar;

for i=1:N_ekf-1
    [x(:,i+1), P_e(:,:,i+1)] = predict(x(:,i), P_e(:,:,i), dt, friction, R, g, l);
    [h, H] = meas_model(x(:,i), dt, Q);
    %Meas mod with theta and theta dot
    %[x(:,i+1), P_e(:,:,i+1)] = update_(x(:,i+1), P_e(:,:,i+1), z_bar(round(i*dt/dt_sim)), Q, h, H);
    %Non-linear measurement model with x and y
    [x(:,i+1), P_e(:,:,i+1)] = update_(x(:,i+1), P_e(:,:,i+1), z_bar(i,:)', Q, h', H);
end     
        

%% iEKF
x_i = zeros(2, N_ekf);
x_i(:,1) = x_start;
%Measurement linear with theta and theta dot
%z_bar = th(:,1) + sqrt(Q).*randn(size(th(:,1)));
%Measurment non-linear as x and y
%z_bar = th(1:dt/dt_sim:end,:) + (sqrt(diag(Q))'.*randn(N_ekf,2));

P_i = zeros(2,2,N_ekf);
P_i(:,:,1) = Sigma_bar;

for i=1:N_ekf-1
    [x_i(:,i+1), P_i(:,:,i+1)] = predict(x_i(:,i), P_i(:,:,i), dt, friction, R, g, l);
    %Linear measurement with theta and theta dot
    %[x_i(:,i+1), P_i(:,:,i+1)] = update_i(x_i(:,i+1), x_i(:,i), P_i(:,:,i+1), z_bar(round(i*dt/dt_sim)), Q, dt);
    %Non-linear measurement model with x and y
    [x_i(:,i+1), P_i(:,:,i+1)] = update_i(x_i(:,i+1), P_i(:,:,i+1), z_bar(i,:)', Q, dt);

end    

%% UKF

x_u = zeros(2, N_ekf);
x_u(:,1) = x_start;
%Measurement linear as theta and theta dot
%z_bar = th(:,1) + sqrt(Q).*randn(size(th(:,1)));
%Measurment non-linear as x and y
%z_bar = th(1:dt/dt_sim:end,:) + (sqrt(diag(Q))'.*randn(N_ekf,2));
%Prediction model
f = @(x) x + [dt*x(2,1); (g/l*sin(x(1,1))*dt - friction*x(2,1))] + sqrt(diag(R)).*ones(2,1);
%Measurement model
%Measrung theta and theta dot => linear measurement model
%h = @(x) x(1,1) + x(2,1)*dt;
%h = @(x) x(1,1);
%Measureing x and y => non-linear measurement model
h = @(x) [l*sin(x(1)), l*cos(x(1))];
Q_u = R;
R_u = Q;
P_u = zeros(2,2,N_ekf);
P_u(:,:,1) = Sigma_bar;

for i=1:N_ekf-1
    %linear meas th th dot
    %z = z_bar(round(i*dt/dt_sim));
    %Non-linear meas
    z = z_bar(i,:);
    [x_u(:,i+1), P_u(:,:,i+1)] = ukf(f, x_u(:,i), P_u(:,:,i), h, z', Q_u, R_u);
end


%% uiKF
if uikf_setting
x_ui = zeros(2, N_ekf);
x_ui(:,1) = x_start;
P_ui = zeros(2,2,N_ekf);
P_ui(:,:,1) = Sigma_bar;

for i=1:N_ekf-1
    %linear meas th th dot
    %z = z_bar(round(i*dt/dt_sim));
    %Non-linear meas
    z = z_bar(i,:);
    [x_ui(:,i+1), P_ui(:,:,i+1)] = uiKF(x_ui(:,i), P_ui(:,:,i), f, z', Q, R, dt);
end
end
%% Error plot
%Linear meas with theta and theta dot
%sim_th = th(1:10:end,:)';
%non-linear meas
sim_th = th_linear(1:dt/dt_sim:end,:)';

err = x - sim_th;
err_i = x_i - sim_th;
err_u = x_u - sim_th;
%err_ui = x_ui - sim_th;

figure('Name', 'Error plot')
subplot(2,3,1)
plot(mod(err(1,:)-pi, 2*pi).^2)
title('Error EKF theta')
subplot(2,3,2)
plot(mod(err_i(1,:)-pi,2*pi).^2)
title('Error iEKF theta')
subplot(2,3,3)
plot(mod(err_u(1,:)-pi,2*pi).^2)
title('Error UKF theta')

subplot(2,3,4)
plot(err(2,:))
title('Error EKF theta dot')
subplot(2,3,5)
plot(err_i(2,:))
title('Error iEKF theta dot')
subplot(2,3,6)
plot(err_u(2,:))
title('Error UKF theta dot')

%% Covariance
P_er = reshape(P_e, [4, N_ekf]);
P_ir = reshape(P_i, [4, N_ekf]);
P_ur = reshape(P_u, [4, N_ekf]);

figure('Name', 'Covariance change')
subplot(2,3,1)
plot(P_er(1,:))
title('EKF cov 1,1')
subplot(2,3,2)
plot(P_ir(1,:))
title('iEKF cov 1,1')
subplot(2,3,3)
plot(P_ur(1,:))
title('UKF cov 1,1')

subplot(2,3,4)
plot(P_er(4,:))
title('EKF cov 2,2')
subplot(2,3,5)
plot(P_ir(4,:))
title('iEKF cov 2,2')
subplot(2,3,6)
plot(P_ur(4,:))
title('UKF cov 2,2')


%% Double plot
%plot_setting = 1;
if plot_setting
dplot = figure('Name', 'Double plot');
for i=1:N_sim
    cla(dplot)
    hold on
    plot(l*sin(x(1,ceil(i/dt*dt_sim))),l*cos(x(1,ceil(i/dt*dt_sim))), '-or','MarkerSize',5,'MarkerFaceColor','r', 'Linewidth', 3)
    plot(l*sin(x_i(1,ceil(i/dt*dt_sim))),l*cos(x_i(1,ceil(i/dt*dt_sim))), '-om','MarkerSize',5,'MarkerFaceColor','m', 'Linewidth', 3)
    plot(l*sin(x_u(1,ceil(i/dt*dt_sim))),l*cos(x_u(1,ceil(i/dt*dt_sim))), '-og','MarkerSize',5,'MarkerFaceColor','g', 'Linewidth', 3)
    
    plot([0, l*sin(th_linear(i,1))],[0, l*cos(th_linear(i,1))], '-ob','MarkerFaceColor','b', 'linewidth', 4)
    hold off
    title(sprintf('Time %f', round(i*dt_sim)))
    axis([-1 1 -1 1])
    pause(dt_sim)
end
end
