%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: To verify the accuracy of the code by inputting a given set of
% intial values and comparing the outputted trajectory and thrust plots
% with those supplied in the "verification case document".
%
% Inputs: None
%
% Outputs: This program calculates the maximum distance and height
% traversed by the bottle rocket using the given set of parameters. It then
% outputs these values to the screen. It also produced two plots, a plot of
% the trajectory of the bottle rocket and a plot of the thrust profile with
% time over the flight, with markers on the plot to indicate where the
% transition between the three phases of flight occur.
%
% Assumptions: It is assumed that the experimental procedure outlined in 
% the report is valid. It is also assumed that ODE45 should be used to
% accomplish the above task. Finally, it is assumed that the verification
% case is accurate.
% 
% Author's ID Number: 60 
% Date Created: 11/26/17
% Date Modified: 12/7/17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% Define global variables
global g % acceleration due to gravity
global gamma % ratio of specific heats
global R % gas constant of air
global m_b % mass of empty 2-liter bottle with cones and fins
global v_b % volume of empty bottle
global c_d % discharge coefficient
global A_t % cross sectional area of throat
global p_w % density of water 
global p_a % ambient air density
global P_a % ambient/atmospheric pressure of air
global Pi_g % initial gage pressure of air in bottle
global Pi_0 % intial total pressure of air in rocket
global Ti_a % initial temperature of the air
global vi_w % initial volume of water inside bottle
global vi_a % initial volume of air inside bottle
global P_end % the pressure of air in the bottle after the water has been expelled
global T_end % the temperature of air in the bottle after the water has been expelled
global mi_a % intial mass of air in the bottle
global C_D % coefficient of drag
global A_B % cross section area of bottle
global x0 % initial x position of the rocket
global z0 % initial height of rocket, ie height of test stand
global l_rod % length of rod on test stand
global theta % initial angle of rocket on test stand
global Vx_0 % intial horizontal velocity of the rocket
global Vy_0 % initial vertical velocity of the rocket
global mi_R % initial mass of rocket

% Set g
g = 9.81; % m/s^2
% Set gamma
gamma = 1.4; % no units
% Set R
R = 287; % J/kgK
% Set m_b
m_b = 0.15; % kg
% Set v_b
v_b = 0.002; % m^3
% Set c_d
c_d = 0.8; % no units
% Set A_t
A_t = pi/4*(2.1/100)^2; % m^2
% Set p_w
p_w = 1000; % kg/m^3
% Set p_a
p_a = 0.961; % kg/m^3
% Set P_a
%P_a = psi2Pa(12.1);
P_a = 83426.56; % Pa
% Set P^i_g
%Pi_g = psi2Pa(50);
Pi_g = 344738; % Pa
% Set P^i_0
Pi_0 = P_a + Pi_g; % Pa
% Set T^i_a
Ti_a = 300; % K
% Set v^i_w
vi_w = 0.001; % m^3
% Set v^i_a
vi_a = v_b-vi_w; % m^3
% Set P_end
P_end = Pi_0*((vi_a/v_b)^gamma); % Pa
% Set T_end
T_end = Ti_a*((vi_a/v_b)^(gamma-1)); % k
% Set mi_a
mi_a = (Pi_0*vi_a)/(R*Ti_a); % kg
% Set C_D
C_D = 0.5;
% Set A_B
A_B = pi/4*(10.5/100)^2; % m^2
% Set x0
x0 = 0; % m
% Set z0
z0 = 0.01; % m 
% Set l_rod
l_rod = 0.5; % m
% Set theta
theta = 45; % degrees
theta = theta * (pi/180); % radians
% Set Vx_0
Vx_0 = 0; % m/s
% Set Vy_0
Vy_0 = 0; % m/s
% Set mi_R
mi_R = m_b+p_w*(v_b-vi_a)+(Pi_0/(R*Ti_a))*vi_a;

% Build initial conditions vector
xi = [x0,z0,Vx_0,Vy_0,vi_a,mi_a,mi_R];

% Call the rocket function to analyze the trajectory of the bottle rocket
[t,x] = ode45('rocket',[0 5],xi);

% Parse output vector for ease of graphing
d = x(:,1); % distance (x position)
h = x(:,2); % height (z position)
Vx = x(:,3); % velocity (x component)
Vz = x(:,4); % velocity (z component)
v_a = x(:,5); % volume of air
m_a = x(:,6); % mass of air
m_R = x(:,7); % mass of rocket

% Call ThrustVec function to produce set of thrust values for plotting
F = ThrustVec(v_a,m_a);

% find the maximum distance
maxd = max(d(h>0))
% find the maximum height
maxh = max(h(h>0))

% Plot Height vs. Distance
figure(1)
plot(d(h>-1),h(h>-1),'LineWidth',2)
title('Height vs Distance');
ylabel('Height (m)');   
xlabel('Distance (m)');

% Plot Thrust vs. Time
figure(2)
hold on 
plot(t(t<=0.5),F(t<=0.5),'LineWidth',2); %t<=0.5 was chosen to best match verification case.
plot(t(find(v_a>v_b,1)),F(find(v_a>v_b,1)),'r*'); % Mark on the graph the transition between phase 1 and phase 2
plot(t(find(F==0,1)),F(find(F==0,1)),'g*'); % Mark on the graph the transition between phase 2 and phase 3
title('Thrust vs Time');
ylabel('Thrust (N)');
xlabel('Time (s)');
legend('Thrust Profile','Transition Between Phase 1 and Phase 2','Transition Between Phase 2 and Phase 3');
hold off