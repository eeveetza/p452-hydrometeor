function R = raincell_structure(rho, Rm)
%aincell_structure Computes rainfallrate as a function of distance from the rain cell center
% as defined in Section 5.3.1.3
%
%     Input parameters:
%     rho        -   Radial distance from the center of the rain cell (km)
%     Rm         -   Peak rain rate at the center  
%
%     Output parameters:
%     R          -   Rainfall rate at distance rho from the center of the rain cell
%
%
%     Example:
%     R = raincell_structure(rho, Rm)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version

if (Rm <=0.4)
    error('Peak rain rate must be greater than 0.4 mm/hr');
end

rho0 = (10.0 - 1.5 * log10(Rm)) / log(Rm / 0.4); %km

% Issue: this equation is valid for Rm > 0.4 mm/hr, but it does not define the
% case Rm <= 0.4

R = Rm * exp(-rho/rho0);

return
end