function [V10, V20, V12] = station_boresight_axes(eps1loc, alpha1loc, eps2loc, alpha2loc, delta)
%station_boresight_axes Computes station antenna boresight axes (main beams)
% as defined in Section 5.2.1.1
%
%     Input parameters:
%     eps1loc    -   Local elevation angle of boresight axis of Station 1 antenna (rad) 
%     alpha1loc  -   Local azimuth angle of boresight axis of Station 1 antenna(rad)
%     eps2loc    -   Local elevation angle of boresight axis of Station 2 antenna (rad) 
%     alpha2loc  -   Local azimuth angle of boresight axis of Station 2 antenna(rad)
%     delta      -   Angle substended by the two stations at the Earth centre (delta = d/ae)
%
%     Output parameters:
%     V10        -   Unit vector for the boresight axis of Station 1 antenna
%     V20        -   Unit vector for the boresight axis of Station 2 antenna
%     V12        -   Unit vector from Station 1 to Station 2
%
%     Example:
%     [V10, V20, V12] = station_boresight_axes(eps1loc, alpha1loc, eps2loc, alpha2loc, delta);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version

% Equation (76)
V10 = [cos(eps1loc).*cos(alpha1loc); ...
       cos(eps1loc).*sin(alpha1loc); ...
       sin(eps1loc)];

% Equation (77)
V20 = [sin(eps2loc).*sin(delta) - cos(eps2loc).*cos(alpha2loc).*cos(delta); ...
       -cos(eps2loc).*sin(alpha2loc); ...
       sin(eps2loc).*cos(delta) + cos(eps2loc).*cos(alpha2loc).*sin(delta)];

% Equation (77b)
V12 = [1; ...
       0; ...
       0];
return
end