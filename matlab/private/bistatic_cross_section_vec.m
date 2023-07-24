function sigma = bistatic_cross_section_vec(eta, phis, alpha_vi, alpha_vs)
%bistatic_cross_section Computes bistatic cross sections
% as defined in Section 5.3.4.4
%
%     Input parameters:
%     eta        -   Raindrop bi-static cross section
%     phis       -   Scattering angle from Station 1 to integration point A
%     alpha_vi   -   The angle that rotates from the incident vertical polarization 
%                    anticlockwise to the polarization perpendicular to the scattering plane
%     alpha_vs   -   The angle that rotates from the scattered vertical polarization 
%                    anticlockwise to the polarization perpendicular to the scattering plane
%
%     Output parameters:
%     sigma      -   bistatic cross section sigma(vv, vh; hv, hh) 2x2 tensor
%
%
%     Example:
%     

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v1    24JUL23     Ivica Stevanovic, OFCOM         Initial version

sigma = zeros(2,2);

sigma(1,1) = eta * (cos(phis)*cos(alpha_vs)*cos(alpha_vi) + sin(alpha_vs)*sin(alpha_vi)).^2;
sigma(1,2) = eta * (cos(phis)*cos(alpha_vs)*sin(alpha_vi) - sin(alpha_vs)*cos(alpha_vi)).^2; 
sigma(2,1) = eta * (cos(phis)*sin(alpha_vs)*cos(alpha_vi) - cos(alpha_vs)*sin(alpha_vi)).^2; 
sigma(2,2) = eta * (cos(phis)*sin(alpha_vs)*sin(alpha_vi) + cos(alpha_vs)*cos(alpha_vi)).^2; 

return
end