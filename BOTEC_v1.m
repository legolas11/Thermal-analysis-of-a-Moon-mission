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

% Steady state nodal calculations
%% Define nodal plan
np.name = {'Shell', 'Structure', 'Payload', 'Engine', 'Battery', 'SP', 'PCB'};
n = length(np.name);
hmatrix = readtable('TMM Spacecraftv2.1.xlsx','Sheet','hmatrix');
np.h = table2array(hmatrix(:,2:end)); %Heat conductance from ith node to jth node, symmetric n*n matrix
np.F = ones(n)./n; %view factor from ith node to jth node, symmetric n*n matrix, sumFij = 1
Nodalplan = readtable('TMM Spacecraftv2.1.xlsx','Sheet','Nodal plan');
np.epsilon = Nodalplan.epsilon; %emittance of nodes
np.epsilonij = ones(n);
for i = 1:n
    for j = 1:n
        np.epsilonij(i,j) = np.epsilon(i)*np.epsilon(j)/(np.epsilon(i)+np.epsilon(j)-np.epsilon(i)*np.epsilon(j));
    end
end
np.alpha = Nodalplan.alpha;
np.A = Nodalplan.A; %surface area of node i
np.Aspace = Nodalplan.Aspace; %effective area with unobstructed view of space
np.Asolar = Nodalplan.Asolar; %effective area with unobstructed view of space
np.Aalbedo = Nodalplan.Aalbedo; %effective area with unobstructed view of space
np.Aplanetary = Nodalplan.Aplanetary; %effective area with unobstructed view of space
np.Qexternal = Nodalplan.Qexternal;
np.Q = Nodalplan.Qinternal_generated_;
np.T0 = Nodalplan.T0; 

%% Set up linear equations system A*T = B
A = zeros(n);
B = zeros(n,1);

for k = 1: 1
    for i = 1:n
        for j = 1:n
            if i == j
                A(i,j) = sum(np.h(i,:)) + 4*sigma*np.T0(i)^3*(np.Aspace(i)*np.epsilon(i)) + sum(np.A(i).*np.F(i,:)'.*np.epsilonij(i,:)')-np.h(i,j)'+4*sigma*np.T0(j).^3.*np.A(i).*np.F(i,j).*np.epsilonij(i,j);
            else
                A(i,j) = -np.h(i,j)'+4*sigma*np.T0(j).^3.*np.A(i).*np.F(i,j).*np.epsilonij(i,j);
            end
    
        end
        B(i,1) = np.Qexternal(i) + np.Q(i) + 3*sigma*np.T0(i)^4*np.Aspace(i)*np.epsilon(i) + 3*sigma*sum((np.T0(i)^4-np.T0(:)'.^4)*np.A(i).*np.F(i,:).*np.epsilonij(i,:));
    
    end

    T = linsolve(A,B);
    
    if abs(sum(T-np.T0)) < 1
        break
    else
        np.T0 = T;
    end
end
