function [x0, y0, h0, d1p, r1p, d2p, r2p] = hydrometeor_geometry_init(d, r1, eps1loc, alpha1loc, h2)
%chydrometeor_geometry_init Constructs hydrometeor scatter geometry
% only the elements that do not depend on the position of the integration element A
% as defined in Section 5.3.3
%
%     Input parameters:
%     d          -   Great circle distance between Stations (km)
%     r1         -   Distance r1 from Station 1 as defined in (80)
%     eps1loc    -   Local elevation angle of boresight axis of Station 1 antenna (rad) 
%     alpha1loc  -   Local azimuth angle of boresight axis of Station 1 antenna(rad)
%     h2         -   Height of the Station 2 above reference plane as defined in Equation (80) (km)
%       
%     Output parameters:
%     x0,y0,h0   -   Center of rain cell (km,km,km)
%     d1p, d2p   -   Horizontal distance from Stations 1,2 to rain cell (km)
%     r1p, r2p   -   Total distance from Stations 1,2 to rain cell (km)
%     Example:
%     [x0, y0, h0, d1p, r1p, d2p, r2p] = hydrometeor_geometry_init(d, r1, eps1loc, alpha1loc, h2);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version


% Flat Earth approximation (Section 5.2.1.1)
eps1 = eps1loc;
alpha1 = alpha1loc;

% Equation (95)
x0 = r1*cos(eps1)*cos(alpha1);
y0 = r1*cos(eps1)*sin(alpha1);
h0 = r1*sin(eps1);

% Distances d1p, d2p, r1p, r2p  from Station 1,2 to the rain cell, Equation (97)

d1p = sqrt(x0^2 + y0^2);
r1p = sqrt(d1p^2 + h0^2);

d2p = sqrt((x0-d)^2 + y0^2);
r2p = sqrt(d2p^2 + (h0-h2)^2);

return
end