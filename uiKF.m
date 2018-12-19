%Merging UKF prediction and iEKF update
function [x, P] = uiKF(x_bar, P_bar, func, z, Q, R, dt)
%UKF prediction step
%Number of states
M_state = numel(x_bar);

%Tunable Paramenters
alpha = 1e-3;
ki = 0;
beta = 2;

lambda = alpha^2*(M_state+ki) - M_state;
c = M_state + lambda;

%Weights
Wm = [lambda/c 0.5/c+zeros(1,2*M_state)];
Wc = Wm;
Wc(1) = Wc(1) + (1-alpha^2+beta);
c = sqrt(c);

%Sigma points
X = sigmas(x_bar, P_bar, c);

%Unscented transform of state
[x1, X1, P1, X2] = ut(func, X, Wm, Wc, M_state, R);

%iEKF update step
[x, P] = update_i(x1, P1, z, Q, dt);

end