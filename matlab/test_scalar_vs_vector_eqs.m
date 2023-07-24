% This script is used to test equations in hydrometeor_geometry
% The vector type equations defined as functions of the angles seem to be
% correct.
% The scalar equations have several issues
% - D_A2 in (97), d_B2 in (101) have a wrong sign in front of the factor cos(alpha2-phi)
% - phi_s in (115) has a wrong sing in the argument of acos()
% - alpha_vi in (124) is incorrect and has several missing factors,
%   introduced a corrected equation based on equation (85) and (86) from the Fascicle 
% - scalar product in equation (126) has a multiplicative factor d missing
% - the same scalar product in (127) is incorrect and not used. Instead (127) is used
% - the scalar product in (129) is incorrect and not used. Instead (128) is used) 

% With the above corrections, both vector and scalar forms give the same results 


f = 12;  %GHz
d = 50;  %km
lat1 = -29.58;  %deg 
lon1 = 30.57;   %deg
lat2 = -29.57;  %deg
lon2 = 30.57;   %deg
h1loc = 8/1000; %km amsl
h2loc = 8/1000; %km amls
alpha1loc = 0    /180.0*pi;  %rad
alpha2loc = 0    /180.0*pi;  %rad
eps1loc = 1      /180.0*pi;   %rad
eps2loc = 38.5   /180.0*pi;   %rad
G1 = 10^(26.3/10);     %linear
G2 = 10^(37.14/10);    %linear
BW1 = 1.5    / 180.0*pi;      %rad       
BW2 = 0.15 / 180.0*pi;      %rad
p = 2; %horizontal
q = 2; %horizontal
press = 1013.25;
wat_vap_dens = 20;
Rm = 10;
T = 20;
M1 = 10;
M2 = 10;
M3 = 10;

% parameter definition

c = 0.23026;

% Calculate the longitude and latitude of the mid-point of the path, phim_e,
% and phim_n for dpnt = 0.5dt
Re = 6371;
dpnt = 0.5*d;
[phim_e, phim_n, ~, ~] = great_circle_path(lon2, lon1, lat2, lat1, Re, dpnt);

% Find radio-refractivity lapse rate DN 
% using the digital maps at phim_e (lon), phim_n (lat) - as a bilinear interpolation

DN = get_interp2('DN50',phim_e,phim_n);

% Compute the effective Earth radius

[ae, ~] = earth_rad_eff(DN);

% The angle subtended by the two stations at the Earth center

delta = d/ae;

% 5.2.1.1 Station antenna boresight axes (main beams)

[V10, V20, V12] = station_boresight_axes(eps1loc, alpha1loc, eps2loc, alpha2loc, delta);

% 5.2.1.3 The off-axis squint angles and related distances

[ksi1, ksi2, rs, r1, r2, h2, relevant] = off_axis_squint_angles(V10, V20, h1loc, h2loc, d, delta);

[x0, y0, h0, d1p, r1p, d2p, r2p] = hydrometeor_geometry_init(d, r1, eps1loc, alpha1loc, h2);

% Pick an arbitrary point in cylindrical coordinates
rho = 0.01;
phi = pi/15;
z = 0.01;


[rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, thetaA2, phis, alpha_vi, alpha_vs] = hydrometeor_geometry(rho, phi, z, d, x0, y0, h0, h2, d1p, r1p, d2p, r2p, ...
               eps1loc, eps2loc, alpha1loc, alpha2loc);


[rA1_v, rA2_v, dB1_v, dB2_v, epsA1_v, epsA2_v, alphaA1_v, alphaA2_v, thetaA1_v, thetaA2_v, phis_v, alpha_vi_v, alpha_vs_v] = hydrometeor_geometry_vec(rho, phi, z, d, x0, y0, h0, h2, eps1loc, eps2loc, alpha1loc, alpha2loc);


fprintf(1,'delta rA1 = %g\n', rA1-rA1_v);
fprintf(1,'delta rA2 = %g\n', rA2-rA2_v);
fprintf(1,'delta dB1 = %g\n', dB1-dB1_v);
fprintf(1,'delta dB2 = %g\n', dB2-dB2_v);
fprintf(1,'delta epsA1 = %g\n', epsA1-epsA1_v);
fprintf(1,'delta epsA2 = %g\n', epsA2-epsA2_v);
fprintf(1,'delta alphaA1 = %g\n', alphaA1-alphaA1_v);
fprintf(1,'delta alphaA2 = %g\n', alphaA2-alphaA2_v);
fprintf(1,'delta thetaA1 = %g\n', thetaA1-thetaA1_v);
fprintf(1,'delta thetaA2 = %g\n', thetaA2-thetaA2_v);
fprintf(1,'delta phis = %g\n', phis-phis_v);
fprintf(1,'delta alpha_vi = %g\n', alpha_vi-alpha_vi_v);
fprintf(1,'delta alpha_vs = %g\n', alpha_vs-alpha_vs_v);


Lpq = tl_p452_hydrometeor(d, f, lat1, lon1, lat2, lon2, h1loc, h2loc, alpha1loc, alpha2loc, eps1loc, eps2loc, G1, G2, BW1, BW2, p, q, press, wat_vap_dens, Rm, T, M1, M2, M3);

Lpq_v = tl_p452_hydrometeor_vec(d, f, lat1, lon1, lat2, lon2, h1loc, h2loc, alpha1loc, alpha2loc, eps1loc, eps2loc, G1, G2, BW1, BW2, p, q, press, wat_vap_dens, Rm, T, M1, M2, M3);


fprintf(1,'delta Lpq = %f\n', Lpq-Lpq_v);