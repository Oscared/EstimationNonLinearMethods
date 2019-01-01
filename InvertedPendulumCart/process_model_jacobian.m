function G = process_model_jacobian(x, dt, constants)
% Function that calculates the numerical jacobian of the process models 
% to play.m
%
% INPUTS:     
%
%   x                (2 or 4)x1     : simulated states
%   dt                    float     : time step
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   G           (2 or 4)x(2 or 4)   : jacobian of process model evaluated 
%                                     at x, size depends on the sim_type 
%                                     and meas_type in play.m


    % step size with the numerical jacobian
    h = 1e-5;
    h_diag = diag(h*ones(size(x)));
    
    G = process_model(repmat(x, size(x')) + h_diag, dt, constants) ...
        - process_model(repmat(x, size(x')) - h_diag, dt, constants);
    G = G./(2*h);

end