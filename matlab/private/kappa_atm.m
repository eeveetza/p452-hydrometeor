function katm = kappa_atm(rA1, rA2, f, p, rho, T)
%kappa_atm Computes atmospheric attenuation
% as defined in Section 5.3.4.2
%
%     Input parameters:$
%     rA1, rA2   -   Distance from Station 1,2 to the integration element A(rho, phi, z)    
%     f          -   Frequency (GHz)
%     p          -   Dry air pressure (hPa)
%     rho        -   Water vapor density (g/m^3)
%     T          -   Temperature (deg C)  
%
%     Output parameters:
%     katm       -   attenuation due to absorption by atmospheric gases
%
%
%     Example:
%     

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version

% issues:
% gamma_atm(t1,2) and gamma_R1,2 (t1,2) are not defined
% they seem to be constant along the path
% in that case equations given in Sections 5.3.4.2 reduce significantly

Tkelvin = 273.15 + T;

[g_0, g_w] = p676d11_ga(f, p, rho, Tkelvin);

gamma_atm = g_0 + g_w; % dB/km

% Equation (118) taking into account that  gamma_atm does not depend on the
% path

katm = gamma_atm * (rA1 + rA2);

return
end