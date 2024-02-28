%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 251 Spring 2024
%
% Assignment Information
%   Assignment:     PM5
%   Authors:        Deepesh Balwani, dbalwani@purdue.edu
%   Team:           R101
%
%   Program Title: Delta V Calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DELTA V1 (Launch from Earth to parking orbit)

% Earth Parameters
earthRotSpeed = 7.292 * 10^(-5);          % Earth's Rotational Speed [rad/s]
earthAvgRadius = 6371;                    % Earth's Average Radius [km]
earthGravConst = 3.986 * 10^5;            % Earth's GM Constant [km^3 / s^2]

% Orbit Parameters
alt = 800;                                % Parking orbit altitude range [km]
orbitRadius = earthAvgRadius + alt;       % Radius of orbit [km]
deltaVLoss = 1.7;                         % DeltaV_Loss [km/s]

% Launch Parameters
latKSC = 28.524;                          % Latitude of GSC [deg]
launchAz = 90;                            % Azimuth angle [deg]

% Calculations (Launch to LEO)
VEH_KSC = earthRotSpeed * earthAvgRadius * cosd(latKSC) * sind(launchAz);
velOrbit = sqrt((earthGravConst) ./ orbitRadius);
deltaVLEO = velOrbit - 0;

% Total DeltaV for this scenario
deltaV1 = deltaVLoss - VEH_KSC + deltaVLEO;


%% DELTA V2 (Plane change of parking orbit to match moon’s orbital plane)

% Plane Change Parameters
earthInc = 23.5;                % Earth's Inclination [deg]
initInc = latKSC - earthInc;    % Initial Orbit Inclination [deg]
finalInc = 5.14;                % Moon's Orbit Relative Inclination [deg]
incDiff = finalInc - initInc;   % Inclination Difference

% Use the V_PC Formula
deltaV2 = 2 * velOrbit * sind(abs(incDiff) / 2);


