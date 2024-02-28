%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 251 Spring 2024
%
% Assignment Information
%   Assignment:     PM5
%   Authors:        Deepesh Balwani, dbalwani@purdue.edu, Jacob Devine, devine38@purdue.edu
%   Team:           R101
%
%   Program Title: Delta V Calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%% DELTA V1 (Launch from Earth to parking orbit)

% Earth Parameters
earthRotSpeed = 7.292 * 10^(-5);          % Earth's Rotational Speed [rad/s]
earthAvgRadius = 6371;                    % Earth's Average Radius [km]
earthGravConst = 3.986 * 10^5;            % Earth's GM Constant [km^3 / s^2]

% Orbit Parameters
%parking orbit should probably go down- can go down as far as 200km perigee
alt = 800;                                % Parking orbit altitude range [km]
earthOrbitRadius = earthAvgRadius + alt;       % Radius of orbit [km]
deltaVLoss = 1.7;                         % DeltaV_Loss [km/s]

% Launch Parameters
latKSC = 28.524;                          % Latitude of GSC [deg]
launchAz = 90;                            % Azimuth angle [deg]

% Calculations (Launch to LEO)
VEH_KSC = earthRotSpeed * earthAvgRadius * cosd(latKSC) * sind(launchAz);
velOrbit = sqrt((earthGravConst) ./ earthOrbitRadius);
deltaVLEO = velOrbit - 0;

% Total DeltaV for this scenario
deltaV1 = deltaVLoss - VEH_KSC + deltaVLEO;


%% DELTA V2 (Plane change of parking orbit to match the Moon's orbital plane)

% Plane Change Parameters
earthInc = 23.5;                % Earth's Inclination [deg]
initInc = latKSC - earthInc;    % Initial Orbit Inclination [deg]
finalInc = 5.14;                % Moon's Orbit Relative Inclination [deg]
incDiff = finalInc - initInc;   % Inclination Difference

% Use the V_PC Formula
deltaV2 = 2 * velOrbit * sind(abs(incDiff) / 2);


%% DELTA V3 (TLI burn, raising apogee to reach the Moon)

%% DELTA V4 (Orbital injection burn to enter lunar orbit)

% Historically, PFS-2 was able to enter a lunar orbit as low as 90km x 130 km and still maintained a stable orbit for 34 days despite the effects of lunar mascons. We will assume a 90 km circular orbit is stable enough for our purposes
moonAvgRadius = 1737.5;
moonGravConst = 4902.8;
perilune = 90;
lunarOrbitRadius = perilune + moonAvgRadius;
lunarOrbitVelocity = sqrt(moonGravConst / lunarOrbitRadius);
%% DELTA V5 (Optional lowering of lunar orbit)

% We will insert directly into our target orbit to minimize dV losses; thus this step is not necessary
% Consider using this option if we end up using a low-thrust transfer stage
deltaV5 = 0;
%% DELTA V6 (Landing delta v estimate)
% Landing spot is at avg radius
seaLevelLandingSemiMajor = lunarOrbitRadius + moonAvgRadius;
maxAltLandingSemiMajor = lunarOrbitRadius + moonAvgRadius + 10.786;
seaLevelLandingVelocity = sqrt(moonGravConst * ((2 / moonAvgRadius) - (1/seaLevelLandingSemiMajor)));
maxAltLandingLandingVelocity = sqrt(moonGravConst * ((2 / (moonAvgRadius + 10.786)) - (1/maxAltLandingSemiMajor)));



