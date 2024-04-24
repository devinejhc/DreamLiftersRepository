%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 251 Spring 2024
%
% Assignment Information
%   Assignment:     PM5
%   Authors:        Deepesh Balwani, dbalwani@purdue.edu
%                   Jacob Devine, devine38@purdue.edu
%                   Stefano Vitello, svitello@purdue.edu
%   Team:           R101
%
%   Program Title: Delta V Calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%% BASE PARAMETERS
% Earth Parameters

earthRotSpeed = 7.292 * 10^(-5);          % Earth's Rotational Speed [rad/s]
earthAvgRadius = 6371;                    % Earth's Average Radius [km]
earthGravConst = 3.986 * 10^5;            % Earth's GM Constant [km^3 / s^2]
g0 = 9.81;                                % Earth's Gravity [m/s^2]

% Moon Parameters

distCOM = 384400;               % Distance between the COMs of the Moon and Earth [km]
moonMass = 7.34767309 * 10^22;  % Mass of the moon [kg]
earthMass = 5.97219 * 10^24;    % Mass of the Earth [kg]
moonVel = 1.023;                % Velocity of the moon relative to Earth [km/s]
moonGravConst = 4.902 * 10^3;   % Gravitational constant of the moon [km^3 / s^2]
moonAvgRadius = 1737.5;         % Average radius of the Moon [km]

%% DELTA V1 (Launch from Earth to parking orbit)

% Orbit Parameters

alt = 275;                                % Parking orbit altitude range [km]
orbitRadius = earthAvgRadius + alt;       % Radius of orbit [km]
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

% Before TLI Orbit Parameters
injAngl = 0;             % Injection angle at perigee [deg]
injVel = 10.88;                 % Injection velocity [km/s]

% TLI Orbit Parameters
lambda = 33.5;              % Arrival angle to the moon [deg]

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

% Calculate deltaV3
deltaV3 = injVel - velOrbit;
%% DELTA V4 (Orbital injection burn to enter lunar orbit)
% Historically, PFS-2 was able to enter a lunar orbit as low as 90km x 130 km and still maintained a stable orbit for 34 days despite the effects of lunar mascons. We will assume a 90 km circular orbit is stable enough for our purposes
lunarOrbitVelocity = sqrt(moonGravConst / r_p);

% Calculate deltaV4
moonOrbitVel = sqrt(moonGravConst / r_p);
deltaV4 = v_p - moonOrbitVel;

%% DELTA V5 (Optional lowering of lunar orbit)
% We will insert directly into our target orbit to minimize dV losses thus periapsis change is not needed
% Consider using this option if we end up using a low-thrust transfer stage
r_p_new = r_p; %New desired perilune radius
deltaV5 = moonOrbitVel - sqrt(2 * moonGravConst / r_p - moonGravConst / ((r_p + r_p_new) / 2)); %Vis Viva for new orbit subtracted from current velocity
%% DELTA V6 (Landing delta v estimate)

% Phase 1 of landing
Phase1Apolune = r_p;
Phase1Perilune = moonAvgRadius;
Phase1SemiMajor = (Phase1Perilune + Phase1Apolune)/2;
Phase1FinalVelocity = sqrt(moonGravConst * ((2/Phase1Apolune) - (1/Phase1SemiMajor)));
Phase1DeltaV = moonOrbitVel - Phase1FinalVelocity;

% Phase 2 of landing
Phase2Velocity = sqrt(moonGravConst * ((2/(Phase1Perilune + .100)) - (1/Phase1SemiMajor)));
Phase2FinalVelocity = .002;
Phase2DeltaV = abs(Phase2FinalVelocity - Phase2Velocity);

% Phase 3 of landing
g = 9.81/6;
burntime = 60;
Phase3DeltaV = burntime * g / 1000;
deltaV6 = Phase1DeltaV + Phase2DeltaV + Phase3DeltaV;

%% TOTAL DELTAV
totalDeltaV = deltaV1 + deltaV2 + deltaV3 + deltaV4 + deltaV5 + deltaV6;
%% Mass Estimation %% —----------------------------------------------------------------------------
%Current ISP and finert values represent nothing, left is last stage, right is launch
Isp = [315.5 465.5]; %Specific Impulse, add more per stage/different ISP
finert = [.23 .1]; %Inert mass fraction of a propulsion system, add more per stage
%RD843 hydrazine for landing, RL10 for kick and capture, Vulcan Centaur for launch and LEO
%Landing and capture stages based on lengthened AVUM, kick stage based on lengthened Centaur
stagecount = length(Isp);
minitial(1:stagecount + 1) = 1000; %Payload mass in kg for the last stage, further generated masses is each subsequent stages payload
dv = [deltaV6 - Phase1DeltaV deltaV5 + deltaV4 + deltaV3 + Phase1DeltaV] * 1000; %Places delta V's into form more usable for loops
LVmaxpay = 27200; %Launch Vehicle max payload (kg)
whilecond = 0;

while whilecond ~= 1
    %Loop generating mass estimates
    if minitial(stagecount + 1) < LVmaxpay %Checks LEO payload mass is less than LV max mass to LEO
        minitial(1) = minitial(1) + 1; %Adds 1kg to the final payload
    else
        whilecond = 1; %Loop conditional is now satisfied
        minitial(1) = minitial(1) - 1; %Subtracts 1kg to the final payload
    end
    for I = 1:1:stagecount %I runs to max number of burns/stages
        mprop(I) = minitial(I) * (exp(dv(I)/(Isp(I) * g0)) - 1) * (1 - finert(I)) / (1-finert(I) * exp(dv(I)/(Isp(I) * g0))); %Estimates propellant mass
        minert(I) = finert(I)/(1-finert(I)) * mprop(I); %Estimates inert mass
        minitial(I+1) =  minert(I) + minitial(I) + mprop(I); %Adds the mass initial of this stage as the payload of the next
    end
    
end
%% Tank sizing
CentaurTankSize = 69.8;  % Centaur tank size in cubic meters
CentaurLength = 12.68 - 1.45; % Centaur tank length
CentaurPropMass = 20830; % Centaur propellant mass in kg
LEMTankSize = 1.906 * 4; % LEM descent stage tank size in cubic meters
LEMPropMass = 8200; % LEM descent stage propellant mass in kg
LEMPressuranceMass = 22;
CentaurLengthPerPropMass = CentaurLength/CentaurPropMass;
TransferStageTankLength = mprop(2) * CentaurLengthPerPropMass;
PressurantPerPropMassHypergol = LEMPressuranceMass/LEMPropMass;
LandingPressurantMass = mprop(1) * PressurantPerPropMassHypergol;
HeliumSpecificVolume = 1/1.71;
PressurantVolume = HeliumSpecificVolume * LandingPressurantMass;
HypergolSpecificVolume = LEMTankSize / LEMPropMass;
LandingStageTankVolume = HypergolSpecificVolume * mprop(1); % in cubic meters

%% TWR calcs
engineThrust = 2452; % Thrust per landing engine, N
descentMass = minitial(2);
thrustRequired = 1.3 * descentMass * g;
engineCount = ceil(thrustRequired / engineThrust);

