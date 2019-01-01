function H = measurement_model_jacobian(x, constants)
% Function with jacobian of measurement models to play.m
%
% INPUTS:     
%
%   x                (2 or 4)x1     : state
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   H        (1,2 or 3)x(2 or 4)    : jacobian of measurement model 
%                                     evaluated at x, size depends on
%                                     the sim_type and meas_type in play.m

    % constants that will be used
    l = constants(2);
    simulation_type = constants(11);
    measurement_type = constants(12);

    % linear measurement (angle + position of cart (if cart is included))
    if measurement_type == 0
        
        % only pendulum
        if simulation_type == 0
        
            H = [1,0];
            
        % with cart
        elseif simulation_type == 1
                
            H = [1, 0, 0, 0;
                 0, 0, 1, 0];
            
        end
    
    % non-linear measurement (x,y for angle + pos of cart (if cart is included))        
    elseif measurement_type == 1    

        % only pendulum
        if simulation_type == 0

            H = [-l*sin(x(1)), 0;
                  l*cos(x(1)), 0];
            
        % with cart
        elseif simulation_type == 1
            
            H = [1, 0, 0, 0;
                 0, 0, -l*sin(x(3)), 0;
                 0, 0,  l*cos(x(3)), 0];
            
        end
    end
end