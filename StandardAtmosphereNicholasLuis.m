%% Functions
function [T,P,rho]=StandardAtmosphereNicholasLuis(h)
g = 9.80667; %acceleration due to gravity (m/s^2)
R = 287.05; %ideal gas const (J/kg*K)

%{
Calculates atmosphere properties for ISA on a standard day.

Inputs
h altitude (in meters)

outputs
T temperature (in K)
P pressure (in Pa)
rho density (in kg/m3)
%}

if and(h>0,h<=11000) % Troposphere
    
    a = -6.5/1000; % lapse rate (K/m)
    T0 = 15.0 + 273.15; % Known sea level temperature (K)
    P0 = 101.325*1000; % Known sea level pressure (Pa)
    rho0 = 1.225; % Known sea level density (kg/m^3)
    h0 = 0; % Altitude of Troposphere base (m)
    
    variableTemp();

    return;

elseif and(h>11000,h<=20000) % we're in the tropopause (const Temp)

    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(11000); % Calculates the values at the top of previous layer (base values of current layer)
    T = T0;
    h0 = 11000;

    constTemp();

    return;

elseif and(h>20000,h<=32000) % we're in the lower stratosphere
    
    a = 1.0/1000; % lapse rate (K/m)
    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(20000); % Calculates the values at the top of previous layer (base values of current layer)
    h0 = 20000;
    
    variableTemp();

    return;

elseif and(h>32000,h<=47000) % we're in the upper stratosphere

    a = 2.8/1000; % lapse rate (K/m)
    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(32000); % Calculates the values at the top of previous layer (base values of current layer)
    h0 = 32000;
    
    variableTemp();
    
    return;

elseif and(h>47000,h<=51000) % we're in the stratopause

    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(47000); % Calculates the values at the top of previous layer (base values of current layer)
    T = T0;
    h0 = 47000;

    constTemp();

    return;

elseif and(h>51000,h<=71000) % we're in the lower mesosphere
% lapse rate is -2.8K/km

    a = -2.8/1000; % lapse rate (K/m)
    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(51000); % Calculates the values at the top of previous layer (base values of current layer)
    h0 = 51000;
    
    variableTemp();
    
    return;

elseif and(h>71000,h<=85000) % we're in the upper mesosphere
% lapse rate is -2.0K/km

    a = -2.0/1000; % lapse rate (K/m)
    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(71000); % Calculates the values at the top of previous layer (base values of current layer)
    h0 = 71000;
    
    variableTemp();
    
    return;

elseif and(h>85000,h<=90000) % we're in the mesopause. The standard atmosphere definition ends at h=86km.

    [T0, P0, rho0] = StandardAtmosphereNicholasLuis(85000); % Calculates the values at the top of previous layer (base values of current layer)
    T = T0;
    h0 = 85000;

    constTemp();

    return;

else % too high or too low! Return not-a-numbers.
    fprintf('Altitude out of bounds!!\n');
    T=NaN;
    P=NaN;
    rho=NaN;
    return;
end

    function constTemp()
        P = P0*exp( -(g/(R*T0))*(h-h0) );
        rho = rho0*exp( -(g/(R*T0))*(h-h0) );
    end

    function variableTemp()
        T = T0 + a*(h-h0);
        P = P0*(T/T0)^(-g/(a*R));
        rho = rho0*(T/T0)^(-g/(a*R) - 1);
    end

    function g = gravity(h)
    %{
    Input: height
    Output: acceleration due to gravity
    %}
        earthRadius = 6378.1*1000; % meters
        g = 9.80665 * (earthRadius / (earthRadius + h))^2;
    
    end

end
