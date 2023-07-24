function krain = kappa_rain(rA1, rA2, f, T, R)
%kappa_atm Computes atmospheric attenuation
% as defined in Section 5.3.4.3
%
%     Input parameters:$
%     rA1, rA2   -   Distance from Station 1,2 to the integration element A(rho, phi, z)    
%     f          -   Frequency (GHz)
%     T          -   Temperature (deg C)
%     R          -   Rainfall rate  
%
%     Output parameters:
%     katm       -   attenuation due to absorption by atmospheric gases
%
%
%     Example:
%     krain = kappa_rain(rA1, rA2, f, T, R);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version

% issues:
% gamma_atm(t1,2) and gamma_R1,2 (t1,2) are not defined
% they seem to be constant along the path
% in that case equations given in Sections 5.3.4.3 reduce significantly

gR = specific_rain_attenuation(f, T, R);

% Equation (118) taking into account that  gamma_atm does not depend on the
% path

krain = gR * (rA1 + rA2);

return
end