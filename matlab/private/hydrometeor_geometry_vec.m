function [rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, thetaA2, phis, alpha_vi, alpha_vs] = hydrometeor_geometry_vec(rho, phi, z, d, x0, y0, h0, h2, eps1loc, eps2loc, alpha1loc, alpha2loc)
%construct_hydrometeor_geometry Constructs hydrometeor scatter geometry
% as defined in Section 5.3.3
%
%     Input parameters:
%     rho,phi,z  -   Cylindrical coordinates (km, rad, km)
%     d          -   Great circle distance between Stations (km)
%     x0,y0,h0   -   Center of rain cell (km,km,km)
%     h2         -   Height of the Station 2 above reference plane as defined in Equation (80) (km)
%     eps1loc    -   Local elevation angle of boresight axis of Station 1 antenna (rad) 
%     alpha1loc  -   Local azimuth angle of boresight axis of Station 1 antenna(rad)
%     eps2loc    -   Local elevation angle of boresight axis of Station 2 antenna (rad) 


%       
%     Output parameters:
%     rA1, rA2   -   Distance from Station 1,2 to the integration element A(rho, phi, z)
%     dB1, dB2   -   Distance from Station 1,2 to the projection of the point A on the ground plane (km)
%     epsA1, epsA2   -   Elevation angles of the position vectors RA1,RA2 (rad)
%     alphaA1, alphaA2   -   Azimuth angles of the position vectors RA1, RA2 (rad)
%     thetaA1, thetaA2   -   Off boresight angles (angles between the unit vectors VA1,2  and the main beam direction of Station 1,2) (rad)   
%     phis               -   Scattering angle from Station 1 to integration point A
%     alpha_vi   -   The angle that rotates from the incident vertical polarization 
%                    anticlockwise to the polarization perpendicular to the scattering plane
%     alpha_vs   -   The angle that rotates from the scattered vertical polarization 
%                    anticlockwise to the polarization perpendicular to the scattering plane

%     Example:
%     [rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, ...
%     thetaA2, phis, alpha_vi, alpha_vs] = hydrometeor_geometry_vec(rho, phi, z, d, x0, y0, h0, h2, eps1loc, eps2loc, alpha1loc);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version


% Flat Earth approximation (Section 5.2.1.1)
eps1 = eps1loc;
eps2 = eps2loc;
alpha1 = alpha1loc;
alpha2 = alpha2loc;

% Equation (98a) 
RA1 = [x0 + rho*cos(phi); y0 + rho*sin(phi); h0 + z];
% Equation (99a) in vector form
rA1 = sqrt(dot(RA1,RA1));
% Equation (64) Fascicle
VA1 = RA1/rA1;

% Equation (65) Fascicle
epsA1 = asin((h0+z)/(rA1));
% Equation (66) Fascicle
alphaA1 = atan2(y0 + rho*sin(phi), x0+rho*cos(phi));

% Equation (68) Fascicle
hA1 =[-sin(alphaA1); cos(alphaA1); 0];
% Equation (70) Fascicle in vector form
dB1 = dot([-(y0 + rho*sin(phi)); x0 + rho*cos(phi); 0], hA1);
% Equation (69) Fascicle
vA1 = [sin(epsA1)*cos(alphaA1); sin(epsA1)*sin(alphaA1); -cos(epsA1)];

% Equation (99a)
RA2 = [x0+rho*cos(phi)-d; y0 + rho*sin(phi); h0 + z - h2];
% Equation (99b) in vector form
rA2 = sqrt(dot(RA2, RA2));
% Equation (74) Fascicle
VA2 = RA2/rA2;

% Before equation (75) Fascicle
% alpha2 = atan2(y0, x0-d);
% Equation (75) Fascicle
epsA2 = asin( (h0-h2+z)/(rA2));
% Equation (76) Fascicle
alphaA2 = atan2(y0 + rho*sin(phi), (x0-d) + rho*cos(phi));

% Equation (78) Fascicle
hA2 = [-sin(alphaA2); cos(alphaA2); 0];

% Equation (81) Fascicle in vector form
dB2 = dot([-(y0 + rho*sin(phi)); x0-d + rho*cos(phi) ; 0], hA2);
%dB2 = sqrt((x0-d + rho*cos(phi))^2 + (y0 + rho*sin(phi))^2);
% Equation (80) Fascicle
vA2 =[sin(epsA2)*cos(alphaA2); sin(epsA2)*sin(alphaA2); -cos(epsA2)];

% Equation (82) Fascicle (vector form)
phis = acos(dot(-VA2, VA1));

% Equation (85) Fascicle
alpha_vi = atan2(dot(vA1, VA2), dot(hA1, VA2));
% Equation (88) Fascicle
alpha_vs = acos(-sin(alpha_vi)*dot(hA1,vA2) + cos(alpha_vi)*dot(vA1,vA2));

% Equations (116) and (117)
thetaA1 = acos(cos(eps1) * cos(epsA1) * cos(alpha1-alphaA1) + sin(eps1)*sin(epsA1));
thetaA2 = acos(cos(eps2) * cos(epsA2) * cos(alpha2-alphaA2) + sin(eps2)*sin(epsA2));



return
end