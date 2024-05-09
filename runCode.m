% Made by Nicholas Luis
clear; clc;

%% Atmosphere Data

% Creating a table of atmosphere values before plotting
dh = 255;
N = round(90000/dh + 1); % Calculates the size of the matrix
AtmosphereValues = NaN(N, 4); % Columns for T, P, rho, and h respectively
n=1;
for h=1:dh:90000
    [Temp, Press, Dens] = StandardAtmosphereNicholasLuis(h);
    AtmosphereValues(n,:) = [Temp, Press, Dens, h];
    n = n+1;
end

% Plotting the Data
figure

subplot(1,3,1)
plot(AtmosphereValues(:,1), AtmosphereValues(:,4))
ylabel("Height (m)");
xlabel("Temperature (K)");

subplot(1,3,2)
plot(AtmosphereValues(:,2), AtmosphereValues(:,4))
xlabel("Pressure (Pa)");

subplot(1,3,3)
plot(AtmosphereValues(:,3), AtmosphereValues(:,4))
xlabel("Density (kg/m^3)");

sgtitle("Atmosphere Values at Various Altitudes");

%% Estimating the Von Karman Line
% Populating the table of values
speedTable = [NaN, NaN, NaN]; % flight speed, orbit speed, and height data

dh = 100;
n=1;
for h=1:dh:90000
    speedTable(n,:) = [FlightSpeed(h), OrbitSpeed(h), h];
    n = n+1;
end

% Plotting the values
figure
plot(speedTable(:,1), speedTable(:,3));
hold on;
plot(speedTable(:,2), speedTable(:,3));
legend("Speed to Maintain Level Flight","Speed to Orbit Earth");
title("Estimating the Von-Karman Line");
ylabel("Height (m)");
xlabel("Speed (m/s)");
hold off;

% FUNCTIONS
function g = gravity(h)
%{
Input: height
Output: acceleration due to gravity
%}
    earthRadius = 6378.1*1000; % meters
    g = 9.8066 * (earthRadius / (earthRadius + h))^2;

end


function orbitV = OrbitSpeed(h)
%{
Input: Height (altitude)
Output: Speed necessary to achieve an orbit

%}
earthRadius = 6378.1*1000; % meters
orbitV = sqrt( gravity(h) * (earthRadius + h)) ;

end


function levelFlightV = FlightSpeed(h)
%{
Input: Height (altitude)
Ouput: Speed necessary to maintain level flight using values of Virgin
       Galactic's SpaceShipOne
%}
m = 1500; % mass of plane (kg)
S = 15; % Wing area (m^2)
Cf = 2; % Lift coefficient

[~, ~, rho] = StandardAtmosphereNicholasLuis(h);

levelFlightV = sqrt( (2*m*gravity(h)) / (Cf*rho*S) );

end
