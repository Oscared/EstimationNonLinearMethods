%Unscented Transform
%Input
%   func = non-linear transform
%   sigma_points = sigma points
%   Wm = weights for mean
%   Wc = weights for covaraince
%   n = number of outputs for func
%   cov = covariance to add
%Output
%   y = transformed mean
%   Y = transformed sigma points
%   cov_bar = transformed covariance
%   diff_bar = transformed differance

function [y, Y, cov_bar, diff_bar] = ut(func, ...
    sigma_points, Wm, Wc, n, cov)
%Number of sigma points
M = size(sigma_points, 2);
Y = zeros(n,M); 
y = zeros(n,1);

for i=1:M
    Y(:,i) = func(sigma_points(:,i));
    y = y + Wm(i)*Y(:,i);
end

diff_bar = Y - repmat(y, [1 M]);
cov_bar = diff_bar*diag(Wc)*diff_bar' + cov;
end