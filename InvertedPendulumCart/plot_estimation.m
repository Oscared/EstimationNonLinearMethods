function plot_estimation(x, color, constants)
% Function for plotting the estimations in play.m
%
% INPUTS:     
%
%   x                (2 or 4)x1     : estimated state
%   color                string     : color specification
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   plots the estimation x in the current figure

    % constants that will be used
    l = constants(2);
    simulation_type = constants(11);
    
    % only pendulum
    if simulation_type == 0
    
        plot(l*sin(x(1)), l*cos(x(1)), 'o', 'Color', color, 'MarkerFaceColor', color, 'linewidth', 4);
        
    
    % with cart    
    elseif simulation_type == 1
        
        plot(x(1) + l*sin(x(3)), l*cos(x(3)), 'o', 'Color', color, 'MarkerFaceColor', color, 'linewidth', 4)

    end
        
end
