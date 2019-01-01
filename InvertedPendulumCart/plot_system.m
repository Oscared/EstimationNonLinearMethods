function plot_system(x, time, constants)
% Function for plotting the pendulum possibly with cart, for play.m
%
% INPUTS:     
%
%   x                (2 or 4)x1     : estimated state
%   time                  float     : current time/time step
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   plots the estimation x in the current figure

    % constatns that will be used
    l = constants(2);
    F = constants(4);
    w = constants(7);
    h = constants(8);
    r = constants(9);
    lim = constants(10);
    simulation_type = constants(11);

    % clear the figure
    h1 = figure(1);
    cla(h1)
    hold on
    
    % only pendulum
    if simulation_type == 0
        
        % pendulum
        plot([0, l*sin(x(1))], [0, l*cos(x(1))],'-ob', 'MarkerFaceColor','b', 'linewidth', 2)
        axis([-l l -l l])
        title(sprintf('Step: %.2f', time))
        grid on
    
    % with cart
    elseif simulation_type == 1    

        % cart
        rectangle('Position',[x(1)-w/2 0-h/2 w h],'Curvature',0.1, 'FaceColor',[0.5, 0.5, 0.5], 'EdgeColor','k')
        rectangle('Position', [x(1)-w/2 0-h r r], 'Curvature',1, 'FaceColor','k','EdgeColor','k')
        rectangle('Position', [x(1)+w/2-r 0-h r r], 'Curvature',1, 'FaceColor','k','EdgeColor','k')

        % pendulum
        plot([x(1), x(1)+l*sin(x(3))], [0, l*cos(x(3))],'-ob', 'MarkerFaceColor','b', 'linewidth', 2)

        % force arrow
        start = x(1) - sign(F)*(w/2 + abs(F)/10);
        stop  = x(1) - sign(F)*w/2;
        plot([start, stop], [0, 0],'r','linewidth', 2);
        plot([stop - sign(F)*0.05, stop], [-0.05, 0],'r','linewidth', 2)
        plot([stop - sign(F)*0.05, stop], [ 0.05, 0],'r', 'linewidth', 2)

        title(sprintf('Step: %.2f', time))
        axis([x(1)-lim-l x(1)+lim+l -lim lim])
        grid on
        
    end 
    
end
