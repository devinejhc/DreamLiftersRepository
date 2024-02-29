%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 251 Spring 2024
%
% Assignment Information
%   Assignment:     PM5
%   Authors:        Deepesh Balwani, dbalwani@purdue.edu
%                   Jacob Devine, devine38@purdue.edu
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
alt = 275;                                % Parking orbit altitude range [km]
orbitRadius = earthAvgRadius + alt;       % Radius of orbit [km]

% parking orbit should probably go down- can go down as far as 200km perigee

deltaVLoss = 1.7;                         % DeltaV_Loss [km/s]

% Launch Parameters
latKSC = 28.524;                          % Latitude of GSC [deg]
launchAz = 90;                            % Azimuth angle [deg]

% Calculations (Launch to LEO)
VEH_KSC = earthRotSpeed * earthAvgRadius * cosd(latKSC) * sind(launchAz);
velOrbit = sqrt((earthGravConst) / orbitRadius);
deltaVLEO = velOrbit - 0;

% Total DeltaV for this scenario
deltaV1 = deltaVLoss - VEH_KSC + deltaVLEO;


%% DELTA V2 (Plane change of parking orbit to match the Moon's orbital plane)

% Plane Change Parameters
earthInc = latKSC;              % Earth's Inclination [deg]
initInc = latKSC - earthInc;    % Initial Orbit Inclination [deg]
finalInc = 5.14;                % Moon's Orbit Relative Inclination [deg]
incDiff = finalInc - initInc;   % Inclination Difference [deg]

% Use the V_PC Formula
deltaV2 = 2 * velOrbit * sind(abs(incDiff) / 2);   


%% DELTA V3 (TLI burn to raise apogee to reach moon)

orbitRadius = 6700;

% Planetary Parameters
distCOM = 384400;               % Distance between the COMs of the Moon and Earth [km]
moonMass = 7.34767309 * 10^22;  % Mass of the moon [kg]
earthMass = 5.97219 * 10^24;    % Mass of the Earth [kg]
moonVel = 1.023;                % Velocity of the moon relative to Earth [km/s]
moonGravConst = 4.90 * 10^3;    % Gravitational constant of the moon [km^3 / s^2]

% Before TLI Orbit Parameters
injAngl = 0;             % Injection angle at perigee [deg]
injVel = 10.88;                 % Injection velocity [km/s]

% TLI Orbit Parameters
lambda = 60;              % Arrival angle to the moon [deg]

% Calculate Sphere of Influence for the moon [km]
radiusInf = distCOM * ((moonMass / earthMass)^(2/5));

% Calculate radius to the Sphere of Influence of the moon [km]
radiusToInf = sqrt(distCOM^2 + radiusInf^2 - 2*distCOM*radiusInf*cosd(lambda));

% Calculate initial specific energy [km^2 / s^2]
initSpecEnergy = (injVel^2 / 2) - (earthGravConst / orbitRadius);

% Calculate initial angular momentum [km^2 / s]
h_0 = orbitRadius .* injVel .* cosd(injAngl); 

% Calculate semi-major axis (a) [km]
a = (-1 * earthGravConst) / (2 * initSpecEnergy);

% Calculate eccentricity of the injection orbit (e) [km]
e = sqrt(1 - (h_0.^2) ./ (earthGravConst * a));

% Calculate arrival phase angle (gamma) [deg]
arrPhaseAngl = asind(radiusInf ./ radiusToInf .* sind(lambda));

% Calculate arrival speed (V1) [km/s]
arrSpeed = sqrt(2 * (initSpecEnergy + earthGravConst / radiusToInf));

% Calculate phi_1 [deg]
phi_1 = acosd(h_0 / (radiusToInf * arrSpeed));

% Calculate periapsis [km]
p = h_0^2 / earthGravConst;

% Calculate theta_1 [deg]
theta_1 = acosd(p/(radiusToInf*e) - (1/e));

% Calculate selenocentric arrival speed (V2) [km/s]
v_2 = sqrt(arrSpeed^2 + moonVel^2 - 2*arrSpeed*moonVel*cosd(phi_1-arrPhaseAngl));

% Calculate alpha and beta of the arrival vector triangle
alpha = acosd((arrSpeed^2 + v_2^2 - moonVel^2) / (2*arrSpeed*v_2));
beta = 180 - alpha - phi_1 + arrPhaseAngl;

% Calculate flight path angle (phi_2) [deg]
phi_2 = 180 - (lambda + beta);

% Calculate arrival trajectory
specEnergy_2 = (v_2^2 / 2) - (moonGravConst / radiusInf);
h_2 = radiusInf * v_2 * cosd(phi_2);
a_2 = (-1 * moonGravConst) / (2 * specEnergy_2);
e_2 = sqrt(1 - (h_2^2)/(moonGravConst*a_2));

% Calculate periselenium radius [km]
r_p = a_2 * (1 - e_2);

% Calculate periapsis speed [km/s]
v_p = sqrt(2 * (specEnergy_2 + (moonGravConst / r_p)));


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



