%% Back of the evelope calculation for thermal analysis for a Moon mission s/c

clear
clc
close all

%% Constants used in calculation

sigma = 5.67*10^(-8);        % Stefan-Boltzmann constant (W/(m^2K^4)
P_sun = 3.856*10^26;         % Total power output from the sun (W)
AU = 149597870700;           % Astronomical unit (m)
R_E = 6378*10^3;             % Radius of Earth

%% Variables choosed manually
a_E = 0.33;                  % Planetary albedo of Earth (0.31 ~ 0.39)
a_M = 0.07;                  % Planetary albedo of Moon 
F_E = 0.5;                   % Visibility factor of Earth, depending on Altitude and angle between local vertical and Sun's ray
alpha = 0.5;                 % Solar absorptance
epsilon = 0.3;               % Infrared emittance
A = 1.5;                     % Surface area of one side of the s/c box-shaped body (m^2)
A_total = 9.5;                 % Total surface area of the s/c

%% Calculation
d = AU;
J_s = P_sun/(4*pi*d^2);      % Solar radiation intensity (W/m^2)
J_a = J_s*a_E*F_E;           % Albedo radiation intensity (W/m^2)
R_orbit_E = R_E + 167000;    % Orbit near Earth
J_p = 237*(R_E/R_orbit_E)^2; % Planetary radiation intensity (W/m^2)
A_s = A;                     % Projected area receiving solar radiation
A_a = A;
A_p = A;
Q = 280;                     % Internally dissipated power (W)

T_E = (A_p*J_p/A_total + Q/(A_total*sigma*epsilon) + (A_s*J_s + A_a*J_a)/(A_total*sigma)*(alpha/epsilon))^(1/4);  %Balanced temperature of the surface near Earth

disp('The balanced temperature (K) of the s/c surface near Earth is:')
disp(T_E)
