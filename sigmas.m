%Sigma points
%Input
%   point = reference point
%   cov = covariance
%   const = coefficient 
%Output
%   sigma_points = sigma points

function [sigma_points] = sigmas(point, cov, const)
change = const*chol(cov)';
fixed_point = point(:,ones(1,numel(point)));
sigma_points = [point fixed_point + change fixed_point - change];
end