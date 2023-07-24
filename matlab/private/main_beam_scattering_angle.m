function phi_ms = main_beam_scattering_angle(V10, V20)
%smain_beam_scattering_angle Computes main beam scattering angle
% as defined in Section 5.2.1.2
%
%     Input parameters:
%     V10        -   Unit vector for the boresight axis of Station 1 antenna
%     V20        -   Unit vector for the boresight axis of Station 2 antenna
%
%     Output parameters:
%     phi_ms     -   Unit vector for the boresight axis of Station 1 antenna
%
%
%     Example:
%     phi_ms = main_beam_scattering_angle(V1, V2);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    20JUL23     Ivica Stevanovic, OFCOM         Initial version

phi_ms = acos(- V20.' * V10);

return
end