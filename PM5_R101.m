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

%% DELTA V1 

% Earth Parameters
earthRotSpeed = 7.292 * 10^(-5);          % Earth's Rotational Speed [rad/s]
earthAvgRadius = 6371;                    % Earth's Average Radius [km]
earthGravConst = 3.986 * 10^5;            % Earth's GM Constant [km^3 / s^2]

% Orbit Parameters
altRange = 800;                           % Parking orbit altitude range [km]
orbitRadius = earthAvgRadius + altRange;  % Radius of orbit [km]
deltaVLoss = 1.7;                         % DeltaV_Loss [km/s]

% Launch Parameters
latKSC = 28.524;                          % Latitude of GSC [deg]
launchAz = 90;                            % Azimuth angle [deg]

% Calculations
VEH_KSC = earthRotSpeed * earthAvgRadius * cosd(latKSC) * sind(launchAz);
velOrbit = sqrt((earthGravConst) ./ orbitRadius);
deltaVLEO = velOrbit - 0;

deltaV1 = deltaVLoss - VEH_KSC + deltaVLEO;
