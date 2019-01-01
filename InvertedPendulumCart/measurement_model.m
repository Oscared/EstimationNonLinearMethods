function h = measurement_model(x, constants)
% Function with measurement models to play.m
%
% INPUTS:     
%
%   x                (2 or 4)x1     : state
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   h                (1,2 or 3)x1   : expected measurement, size depends on
%                                     the sim_type and meas_type in play.m

    % constants that will be used
    l = constants(2);
    simulation_type = constants(11);
    measurement_type = constants(12);

    % linear measurement (angle + position of cart (if cart is included))
    if measurement_type == 0
        
            h = x(1:2:end);
    
    % non-linear measurement (x,y for angle + pos of cart (if cart is included))        
    elseif measurement_type == 1
        
        % only pendulum
        if simulation_type == 0
            
            h = [l*cos(x(1)); 
                 l*sin(x(1))];
        
        % with cart     
        elseif simulation_type == 1
           
            h = [x(1);
                 x(1) + l*cos(x(3));
                 x(1) + l*sin(x(3))];
            
        end     
    
    end     

end