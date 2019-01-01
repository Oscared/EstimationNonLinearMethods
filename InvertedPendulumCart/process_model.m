function process = process_model(x, dt, constants)
% Function with process models to play.m
%
% INPUTS:     
%
%   x                (2 or 4)x1     : simulated states
%   dt                    float     : time step
%   constants              1x12     : useful constants from play.m
%
% OUTPUTS:
%
%   process            (2 or 4)x1   : predicted new state, size depends on
%                                     the sim_type and meas_type in play.m

    % constants that will be used
    g = constants(1);
    l = constants(2);
    f = constants(3);
    F = constants(4);
    m = constants(5);
    M = constants(6);
    simulation_type = constants(11);
    
    % only pendulum
    if simulation_type == 0
        
        process = x + [x(2,:);
                (g/l)*sin(x(1,:))]*dt;
        
        % lower the angular velocity due to friction
        process(2,:) = process(2,:) - x(2,:)*f;
    
    % with cart
    elseif simulation_type == 1
    
        process = x + [x(2,:);
        
                  (-m*g*sin(x(3,:)).*cos(x(3,:)) + m*l*x(4,:).^2.*sin(x(3,:)) + ...
                   f*m*x(4,:).*cos(x(3,:)) + F)./(M + (1-cos(x(3,:)).^2)*m);
               
                   x(4,:);
                   
                   ((M + m)*(g*sin(x(3,:)) - f*x(4,:)) - (l*m*x(4,:).^2.*sin(x(3,:)) + F).* ...
                   cos(x(3,:)))./(l*(M + (1 - cos(x(3,:)).^2)*m))]*dt;  
               
    end
               
end
