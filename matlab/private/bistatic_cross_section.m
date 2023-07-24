function sigma = bistatic_cross_section(eta, phis, alpha_vi, alpha_vs)
%bistatic_cross_section Computes bistatic cross sections
% as defined in Section 5.3.4.4
%
%     Input parameters:
%     rho,phi,z  -   Cylindrical coordinates (km, rad, km)
%     eta        -   Raindrop bi-static cross section
%     phis       -   Scattering angle from Station 1 to integration point A
%     d          -   Great circle distance between Stations 1 and 2 (km)
%     x0,y0,h0   -   Cartesian coordinates of the integration element A (km,km,km)
%     h2         -   Height of the Station 2 above reference plane as defined in Equation (80) (km)
%     rA1, rA2   -   Distance from Station 1,2 to the integration element A(rho, phi, z)
%     dB1, dB2   -   Distance from Station 1,2 to the projection of the point A on the ground plane (km)

%
%     Output parameters:
%     sigma      -   bistatic cross section sigma(vv, vh; hv, hh) 2x2 tensor
%
%
%     Example:
%     

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v1    20JUL23     Ivica Stevanovic, OFCOM         Initial version


% Equations 123


sigma = zeros(2,2);

sigma(1,1) = eta * (cos(phis)*cos(alpha_vs)*cos(alpha_vi) + sin(alpha_vs)*sin(alpha_vi)).^2;
sigma(1,2) = eta * (cos(phis)*cos(alpha_vs)*sin(alpha_vi) - sin(alpha_vs)*cos(alpha_vi)).^2; 
sigma(2,1) = eta * (cos(phis)*sin(alpha_vs)*cos(alpha_vi) - cos(alpha_vs)*sin(alpha_vi)).^2; 
sigma(2,2) = eta * (cos(phis)*sin(alpha_vs)*sin(alpha_vi) + cos(alpha_vs)*cos(alpha_vi)).^2; 

return
end