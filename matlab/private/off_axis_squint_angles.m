function [ksi1, ksi2, rs, r1, r2, h2, relevant] = off_axis_squint_angles(V10, V20, h1loc, h2loc, d, delta)
%off_axis_squint_angles Computes the off axis squint angles
% as defined in Section 5.2.1.3
%
%     Input parameters:
%     V10        -   Unit vector for the boresight axis of Station 1 antenna
%     V20        -   Unit vector for the boresight axis of Station 2 antenna    
%     h1loc      -   Local height above mean sea level of Station 1 (km)
%     h2loc      -   Local height above mean sea level of Station 2 (km)
%     d          -   Great circle distance between Statio n1 and Station 2 (km)
%     delta      -   Angle subtended by the two stations at the Earth centre 
%
%     Output parameters:
%     ksi1,2     -   Off-Axis squint angles at Stations 1 and 2
%     rs,r1,r2   -   Distances (km)
%     h2         -   Height of station 2 above the reference plane
%     relevant   -   boolean flag set to false if phi_ms < 0.001 rad
%
%
%     Example:
%     [ksi1, ksi2] = off_axis_squint_angles(V10, V20, h1loc, h2loc, d, delta);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version

relevant = true;

% Equation 80 

h2 = h2loc - h1loc - d * delta/2;    % km

% Equation 78
% if phi_ms < 0.001 rad, the two antenna beams are eigher approximately
% parallel or co-linear, in which case coupling by rain scatter is
% negligible and there is no need to calculate

phi_ms = main_beam_scattering_angle(V10, V20);

if (phi_ms < 0.001)
    relevant = false;
    return
end

% Equation 81

Vso = cross(V20, V10)./sin(phi_ms);

% Equation 80

R = [Vso, V10, -V20] \ [d; 0; h2];

% Is it possible for rs, r1 and r2 to be negative and how to prevent that
% from happening?
rs = (R(1));
r1 = (R(2));
r2 = (R(3));

% Equation (79)
% if squint angles are less than 3dB BW of the relevant antenna, main beam
% to main beam coupling is possible and calculating hydrometeor scatter is
% required
ksi1 = atan2(rs,r1);
ksi2 = atan2(rs,r2);


return
end