function z = zeta(h, hR)
%zeta Height dependence of the radar reflectivity
%     Input parameters:
%     h       -   Height (km)
%     hR      -   Rain height (km)
%
%     Output parameters:
%     z       -   Height dependence of the radar reflectivity according to ITU-R P.452-18
%
%     Example:
%     z = zeta(h, hR);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version

z = 1;

if (h > hR)
    z = 10.^(-0.65.*(h-hR));
end

return
end