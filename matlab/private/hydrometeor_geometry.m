function [rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, thetaA2, phis, alpha_vi, alpha_vs] = hydrometeor_geometry(rho, phi, z, d, x0, y0, h0, h2, d1p, r1p, d2p, r2p, eps1loc, eps2loc, alpha1loc, alpha2loc)
%construct_hydrometeor_geometry Constructs hydrometeor scatter geometry
%
%     Input parameters:
%     rho,phi,z  -   Cylindrical coordinates (km, rad, km)
%     d          -   Great circle distance between Stations (km)
%     x0,y0,h0   -   Center of rain cell (km,km,km)
%     d1p, d2p   -   Horizontal distance from Stations 1,2 to rain cell (km)
%     r1p, r2p   -   Total distance from Stations 1,2 to rain cell (km)
%     eps1loc    -   Local elevation angle of boresight axis of Station 1 antenna (rad) 
%     alpha1loc  -   Local azimuth angle of boresight axis of Station 1 antenna(rad)
%     eps2loc    -   Local elevation angle of boresight axis of Station 2 antenna (rad) 
%     alpha2loc  -   Local azimuth angle of boresight axis of Station 2 antenna(rad)
%     h2         -   Height of the Station 2 above reference plane as defined in Equation (80) (km)
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

%
%     Example:
%     [rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, thetaA2, phis] = hydrometeor_geometry(rho, phi, z, d, r1, eps1loc, eps2loc, alpha1loc, alpha2loc);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version


% Flat Earth approximation (Section 5.2.1.1)
eps1 = eps1loc;
eps2 = eps2loc;

alpha1 = atan2(y0, x0);
alpha2 = atan2(y0, x0-d);

% Distances rA1, rA2 from Station 1,2 to the integration element A(rho,phi,z), Equation (97)

DA1 = sqrt( 1 + (rho^2 + z^2 + 2*rho*d1p*cos(alpha1-phi)+ 2*h0*z) / (r1p^2) );
rA1 = r1p*DA1;

% ISSUE: Corrected the sign in front of the term 2*rho*d2p
DA2 = sqrt(1 + (rho.^2 + z.^2 + 2*rho*d2p*cos(alpha2-phi) + 2*(h0-h2)*z) /(r2p.^2) );
rA2 = r2p*DA2;

% Equation (101)

dB1 = sqrt(d1p^2 + rho^2 + 2*rho*d1p*cos(alpha1-phi));
dB2 = sqrt(d2p^2 + rho^2 + 2*rho*d2p*cos(alpha2-phi));

% Elevation angles of the position vectors RA1,RA2, Equation 100

epsA1 = atan2((h0+z),dB1);
epsA2 = atan2((h0-h2+z),dB2);

% the azimuth angles of the position vectors RA1, RA2, Equation 102

alphaA1 = atan2( (y0 + rho*sin(phi)) , (   x0 + rho*cos(phi)) );
alphaA2 = atan2( (y0 + rho*sin(phi))  , (x0-d + rho*cos(phi)) ); 

% The scattering angle phis from Station 1 to the integration point A,
% Equation (115)

phis = acos( -rA1/rA2 + ( d* (x0+rho*cos(phi)) + h2*(h0+z) )/(rA1*rA2) );

% Equation (128)
scalar_nu_nu = ( dB1.^2 * rA2.^2 + (h0 - h2 + z) * ( h2*dB1.^2 - (h0+z)*(x0+rho*cos(phi))*d ))/ ...
            ( rA1*dB1*rA2*dB2 );

% Equation (126)
scalar_h_nu = ( d* (h0 - h2 + z)* (y0 + rho * sin(phi) ) ) / (dB1*dB2*rA2);

% Equation 124
x = x0 + rho*cos(phi);
y = y0 + rho*sin(phi);
h = h0 + z;

alpha_vi = atan2( ( h * x * (x-d) + h*y.^2 -dB1^2*(h-h2) ) , (rA1*d*y) );

% Equation 125

alpha_vs = acos( -sin(alpha_vi) * scalar_h_nu + cos(alpha_vi) * scalar_nu_nu  );

% Off boresight angles (angles between the unit vectors VA1,2  and the main
% beam direction of Station 1,2, Equations 116 and 117

thetaA1 = acos( cos(eps1)*cos(epsA1)*cos(alpha1-alphaA1) + sin(eps1)*sin(epsA1) );
thetaA2 = acos( cos(eps2)*cos(epsA2)*cos(alpha2-alphaA2) + sin(eps2)*sin(epsA2) );


return
end