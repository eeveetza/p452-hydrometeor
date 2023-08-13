function Lpq = tl_p452_hydrometeor(d, f, lat1, lon1, lat2, lon2, h1loc, h2loc, alpha1loc, alpha2loc, eps1loc, eps2loc, G1, G2, BW1, BW2, p, q, press, wat_vap_dens, Rm, T, M1, M2, M3)
%tl_p452_hydrometeor transmission loss due to hydrometeor scatterrer according to ITU-R P.452-18
%   Lpq = tl_p452_hydrometeor
%
%   This is the MAIN function that computes the transmission due to hydrometeor scatterer
%   as defined in Rec. ITU-R P.452-18 (Section 5).
%   Other functions called from this function are in ./private/ subfolder.
%
%     Input parameters:
%     d             -   Great-circle path distance (km)
%     f             -   Frequency  (GHz)
%     lat1, lon1    -   Latitude and Longitude for Station 1 (deg)
%     lat2, lon2    -   Latitude and Longitude for Station 2 (deg)
%     h1loc         -   Local heights above mean sea level of Stations 1 and 2 (km)
%     h2loc  
%     alpha1loc     -   Local bearings (azimuth) of Station 1 from Station 2 and
%     alpha2loc         Station 2 from Station 1, in the clockwise sense (rad) 
%     eps1loc       -   Local horizon angles (elevation) for Station 1 and 2 (rad) 
%     eps2loc
%     G1, G2        -   Gain for each antenna as a function of both antenna
%                       boresight angle and antenna polarization (linear!)
%     BW1, BW2      -   Antenna beam widths (either main beam or side lobes)
%                       depending on the required coupling (rad)
%     p, q          -   Polarization of Station 1 (p) and Station 2 (q)
%                       antenna (1 = vertical, 2 = horizontal)
%     press         -   Surface pressure (default 1013.25 hPa)
%     wat_vap_dens  -   Surface water-vapour density (default 8g/m^3)
%     Rm            -   Peak rain rate at the centre of the rain cell (> 0.4 mm/hr)
%     T             -   Surface temperature (default 15 deg C)
%     M1,M2,M3      -   Number of integration points in the numerical
%                       integration over rho, phi, and z, respectively


%
%     Output parameters:
%     Lpq           -   transmission loss due to hydrometeor scatterer according to Rec. ITU-R P.452-18
%
%     Example:
%     Lpq = tl_p452_hydrometeor(d, f, lat1, lon1, lat2, lon2, h1loc, h2loc,alpha1loc, alpha2loc, eps1loc, eps2loc, G1, G2, BW1, BW2, p, q, press, wat_vap_dens, Rm, T, M1, M2, M3);


%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version

% MATLAB Version '9.12.0.1975300 (R2022a) Update 3' used in development of this code
%
% The Software is provided "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED,
% INCLUDING BUT NOT LIMITED TO, THE WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE AND NON-INFRINGMENT OF INTELLECTUAL PROPERTY RIGHTS
%
% Neither the Software Copyright Holder (or its affiliates) nor the ITU
% shall be held liable in any event for any damages whatsoever
% (including, without limitation, damages for loss of profits, business
% interruption, loss of information, or any other pecuniary loss)
% arising out of or related to the use of or inability to use the Software.
%
% THE AUTHOR(S) AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
%
% This function calls other functions that are placed in the ./private folder


% TODO: verify input argument values and limits

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

if(~relevant)
    warning('The rain scatter is negligible. It is not calculated.')
    Lpq = 0;
    return
end

if (ksi1 > BW1 || ksi2 > BW2)
    warning('Off-axis squint angles are larger than the relevant antenna beamwidths. The hydrometeor scatter is negligible and it is not calculated.')
    Lpq = 0;
    return
end

% 5.3.3 Initial hydrometeor geometry, independent of integration, Equation 97

[x0, y0, h0, d1p, r1p, d2p, r2p] = hydrometeor_geometry_init(d, r1, eps1loc, alpha1loc, h2);


% 5.3.1.4 Compute the rain height (amsl)
% ISSUE: for which lat/lon? not defined in Recommendation
% using the mid point
lat = phim_n;
lon = phim_e;

hR = rain_height(lat, lon);

% 5.3.3 Step 3 Compute the integration limits

[rhomax, zmin, zmax] = integration_limits(d, r1, hR, eps1loc, eps2loc, alpha1loc, BW1, BW2);


% 5.3.5 Step 5: Integration of the scatter transfer function

% 1) Determine the Gauss quadrature points for ksi, eta and zeta in [-1,1]
% ISSUE: Step 5.1) says that ksi_i and eta_j need to be computed within the
% loop for n (z_n), but is this necessary? z_n is not dependent on rho_max.
% Also there is a step misssing to determine % the Gauss quadrature for zeta_n, w_n.
% Finally, paragraph 5.3.7 in Step 1) seems to be misreferenced.

[ksi_i,  w_i] = gaussian_quadrature(M1);
[eta_j,  w_j] = gaussian_quadrature(M2);
[zeta_n, w_n] = gaussian_quadrature(M3);

% 2) Introduce the Guass quadrature points into (131) to calculate the
% radii rho_i, the azimuth angles phi_j, and the heights z_n within the
% integration volume

rho_i = rhomax*(ksi_i + 1)/2.0;
phi_j = pi*(eta_j+1.0);
z_n = (zmax-zmin)/2.0 * zeta_n + (zmax + zmin)/2.0;


H_n = zeros(M3,1);

for n = 1 : M3
    % 3) Start with the lower layer within the rain cell n = 1 with the
    % quadrature node zeta_n(1) and weight w_n(1)

    for i = 1 : M1
        for j = 1 : M2

            % 4) Use the resultant radii, azimuth angles and height (rho_i, phi_j,
            % zeta_n) to calculate values of the parameters reported in equations
            % (97)-(102)

            % 5) For each point, use the above values to determine the
            % off boresight angles (116) and (117)

            [rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, thetaA2, phis, alpha_vi, alpha_vs] = hydrometeor_geometry(rho_i(i), phi_j(j), z_n(n), d, x0, y0, h0, h2, d1p, r1p, d2p, r2p, ...
               eps1loc, eps2loc, alpha1loc, alpha2loc);
            % old = phis;
            %[rA1, rA2, dB1, dB2, epsA1, epsA2, alphaA1, alphaA2, thetaA1, thetaA2, phis, alpha_vi, alpha_vs] = hydrometeor_geometry_vec(rho_i(i), phi_j(j), z_n(n), d, x0, y0, h0, h2, eps1loc, eps2loc, alpha1loc);
            %old-phis
            % and hence the (linear!) gain of each antenna
            %TODO: include antenna gains as functions G1,2(thetaA1,2)

            % the atmospheric attenuation (118), and rain attenuation (119).
            
            katm = kappa_atm(rA1, rA2, f, press, wat_vap_dens, T);
            
            R = raincell_structure(rho_i(i), Rm);

            krain = kappa_rain(rA1, rA2, f, T, R);
            
            kappa = katm + krain;

            % the bi-static cross sections, equations (123a)-(123d)

            eta = raindrop_bscross_section(phis, f, R, T);

            sigma = bistatic_cross_section(eta, phis, alpha_vi, alpha_vs);
            %sigma = bistatic_cross_section_vec(eta, phis, alpha_vi, alpha_vs);

            % 6) Use the outcomes of procedure 5 in calculating the
            % corresponding function Fijn(ksi_i, eta_j, zeta_n) given in (104)

            % Height dependence of the radar reflectivity
            h = h0 + z_n(n);  % Equation (96)
            zz = zeta(h, hR); % Equation (75)

            % Equation (104)
            Fijn = G1 * G2 * sigma(p,q) * exp(-c*kappa) * zz / (rA1.^2 * rA2.^2);
          

            % 7) Multiply each Fijn by the corresponding Gauss weights
            % ISSUE: Step 7 suggests to use w_n(1) here, but we are the
            % multiplication with w_n(n) in Step 8 again!
            % Fijn = w_i(i) * w_j(j) * w_n(n) * Fijn;
            Fijn = w_i(i) * w_j(j) * Fijn;

            % 8) Sum all values of Fijn
            H_n(n) = H_n(n) + Fijn;
        end
    end
    % 8) and multiply the result by w_n(n)*pi*(zmax-zmin)*dc^2 yielding Hn
    % ISSUE: dc is not defined, supposing this is the great circle distance d
    H_n(n) = H_n(n) * w_n(n)*pi*(zmax-zmin) * d^2;
    
    % 9) Repeat procedures 1 through 8 with increasing order of n (1:M3)
    % to obtain all the values of Hn as given in (133)
    % ISSUE it should be procedures 3 to 8

end % loop over n

% 10) Sum all values of H_n and divide the the resultant over (r1*r2')^2
% to get the value of the scatter transfer function Cpq as given in (132)

Cpq = sum(H_n)/(r1*r2p).^2;

% 5.1 Transmission loss due to hydrometeor scatterer, Equation 73

Lpq = 73.4399 + 20*log10(f) - 10*log10(Cpq);



return
end