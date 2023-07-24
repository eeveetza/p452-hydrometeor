function [rhomax, zmin, zmax] = integration_limits(d, r1, hR, eps1loc, eps2loc, alpha1loc, BW1, BW2)
%integration_limits Computes integration limits for  hydrometeor scatter geometry
% as defined in Section 5.3.3
%
%     Input parameters:
%     d          -   Great circle distance between Stations (km)
%     r1         -   Distance r1 from Station 1 as defined in (80)
%     hR         -   Rain height (km) as defined in (90)
%     eps1loc    -   Local elevation angle of boresight axis of Station 1 antenna (rad) 
%     eps2loc    -   Local elevation angle of boresight axis of Station 2 antenna (rad) 
%     alpha1loc  -   Local azimuth angle of boresight axis of Station 1 antenna (rad)
%     BW1, BW2   -   Antenna beamwidths of Station 1,2
%
%     Output parameters:
%     rhomax     -   Radial limit (km)
%     zmin, zmax   -   z-axis limits (km)

%     Example:
%     [rhomax, zmin, zmax] = integration_limits(d, r1, eps1loc, eps2loc, BW1, BW2)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version

% Issues in equations (106)-(114)
% eps12 -> vareps1,2
% Equations 108 and 109 should have BW2

% Flat Earth approximation (Section 5.2.1.1)
eps1 = eps1loc;
eps2 = eps2loc;
alpha1 = alpha1loc;


% Equations (95)
x0 = r1*cos(eps1)*cos(alpha1);
y0 = r1*cos(eps1)*sin(alpha1);
h0 = r1*sin(eps1);

% Equations 106-111

z1max = sqrt(x0^2+y0^2) * tan(eps1 + 0.5*BW1)-h0;
z1min = sqrt(x0^2+y0^2) * tan(eps1 - 0.5*BW1)-h0;

z2max = sqrt((x0-d)^2+y0^2) * tan(eps2 + 0.5*BW2)-h0;
z2min = sqrt((x0-d)^2+y0^2) * tan(eps2 - 0.5*BW2)-h0;
zmax = min(max(z1max,z2max), hR);
zmin = min(z1min,z2min);

% Equations 112-114

rho1 =     x0 - d/(1 + tan(eps1 + 0.5*BW1)*cot(eps2-0.5*BW2));
rho2 = d - x0 - d/(1 + tan(eps2 + 0.5*BW2)*cot(eps1-0.5*BW1));

rhomax = 0.5*(rho1+rho2);

return
end