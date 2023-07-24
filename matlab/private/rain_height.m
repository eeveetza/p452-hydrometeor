function hR = rain_height(lat, lon)
%rain_height Computes rain height (amsl) 
% as defined in Section 5.3.1.4
%
%     Input parameters:
%     lat        -   Latitude (deg)
%     lon        -   Longitude (deg) 
%
%     Output parameters:
%     hR         -   rain height (km) amsl
%
%
%     Example:
%     hR = rain_height(lat, lon);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version

hiso = get_interp2('h0',lon, lat);

hR = hiso + 0.36; 

return
end