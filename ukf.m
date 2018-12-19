%UKF
%Input
%   fstate = function of state (prediction)
%   x = state estimate
%   P = state covariance
%   hmean = function of measurement
%   z = measurement
%   Q = process noise covaraince
%   R = measurement noise covariance
%Output
%   x = estimated state
%   P = estimated state covariance

function [x, P] = ukf(fstate, x, P, hmeas, z, Q, R)
%Number of states
M_state = numel(x);
%Number of measurements
m_meas = numel(z);

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
X = sigmas(x, P, c);

%Unscented transform of state
[x1, X1, P1, X2] = ut(fstate, X, Wm, Wc, M_state, Q);
X1 = sigmas(x1, P1, c);
X2 = X1 - repmat(x1, [1 size(X1,2)]); 
%Unscented transform of measurement
[z1, Z1, P2, Z2] = ut(hmeas, X1, Wm, Wc, m_meas, R);

%State-measurement cross-covariance matrix
P12 = X2*diag(Wc)*Z2';

%Kalman gain
K = P12/P2;

%State update
x = x1 + K*(z - z1);

%Covariance update
P = P1 - K*P12';
end
