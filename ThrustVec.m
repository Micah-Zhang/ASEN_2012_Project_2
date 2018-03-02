%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: To calculate the thrust produced by a rocket given a set of
% values for the volume of air inside the bottle and the mass 
% of the air inside the bottle. These calculated thrust values will then be
% used to plot the thrust of the rocket vs. time by the the main
% script.
%
% Inputs: This program takes two inputs, the volume of the air inside the 
% rocket at a given time and the mass of air inside the rocket.
% 
% Outputs: This function uses the numerous equations provided in the project
% description along with inputted values for the volume of air inside the
% rocket and the mass of air inside the rocket. It then uses these values 
% to determine which phase of flight the rocket is in and calculate the 
% thrust produced by the rocket accordingly. This calculated thrust is
% given as output to be used by the main script to create a thrust profile
% for the rocket.
%
% Assumptions: It is assumed that the experimental procedure outlined in 
% the report is valid.
% 
% Author's ID Number: 60 
% Date Created: 11/26/17
% Date Modified: 12/7/17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F] = ThrustVec(vair,mair)
% Import all necessary global variables
global gamma % ratio of specific heats
global R % gas constant of air
global v_b % volume of empty bottle
global c_d % discharge coefficient
global A_t % cross sectional area of throat
global P_a % ambient/atmospheric pressure of air
global Pi_0 % intial total pressure of air in rocket
global vi_a % initial volume of air inside bottle
global P_end % the pressure of air in the bottle after the water has been expelled
global mi_a % intial mass of air in the bottle

for i = 1:length(vair)
    % Calculate pressure of air inside rocket assuming thrust stage 2 or
    % ballistic phase. This pressure will be compared to P_a to determine
    % whether in thrust stage 2 or ballistic phase.
    P = ((mair(i)/mi_a)^gamma)*P_end;
    
    % Determine if in thrust phase 1
    if vair(i) < v_b
        % Calculate the pressure inside the bottle
        P = ((vi_a/vair(i))^gamma)*Pi_0;
        % Find the thrust produced by the rocket in phase 1
        F(i) = 2*c_d*(P-P_a)*A_t;
    % Determine if in thrust phase 2
    elseif P > P_a
        % Find density of air inside the bottle
        p = mair(i)/v_b;
        % Find temperature of air inside the bottle
        T = P/(p*R);
        % Find critical pressure
        P_star = P*((2/(gamma+1))^(gamma/(gamma-1)));
        % Determine if flow is choked
        if P_star > P_a
            % Find the exit pressure
            P_e = P_star;
            % Find the exit temperature
            T_e = (2/(gamma+1))*T;
            % Find the exit density
            p_e = P_e/(R*T_e);
            % Find exit velocity
            V_e = sqrt(gamma*R*T_e);
        else
            % The flow is not choked
            P_e = P_a;
            % Find the exit Mach number
            M_e = sqrt((nthroot(P/P_a,gamma/(gamma-1))-1)/((gamma-1)/2));
            % Find the exit temperature
            T_e = T*(1+((gamma-1)/2)*M_e^2);
            % Find the exit density
            p_e = P_a/(R*T_e);
            % Find the exit velocity
            V_e = M_e*sqrt(gamma*R*T_e);
        end
        % Find dmairdt, the rate of change of the mass of air in the rocket
        mair_dot = c_d*p_e*A_t*V_e;
        % Find the thrust of the bottle rocket
        F(i) = mair_dot*V_e+(P_e-P_a)*A_t;
    else
        % Rocket is in ballistic phase
        F(i) = 0; % In the ballistic phase, no thrust is being produced. Hence the name.
    end    
end
end

