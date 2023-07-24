function gR = specific_rain_attenuation(f, T, R)
%specific_rain_attenuation Computes specific rain attenuation
% as defined in Section 5.3.1.2
%
%     Input parameters:
%     f          -   Frequency (GHz)   x
%     T          -   Temperature (deg C)
%     R          -   Rainfall rate  
%
%     Output parameters:
%     gR         -   Specific rain attenuation (dB/km)
%
%
%     Example:
%     gR = specific_rain_attenuation(f, T, R);

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    21JUL23     Ivica Stevanovic, OFCOM         Initial version

% Table 6

cmi = [0.86481	0.0025984	-3.2727e-05
-0.32507	-0.025593	0.00040852
0.70075	0.041632	-0.00084479
-0.4162	-0.023144	0.00063446
0.11971	0.0054147	-0.00022071
-0.018495	-0.00049312	3.6339e-05
0.0012143	8.1571e-06	-2.2949e-06];

% Table 7

dmi = [-9.2859	-0.026677	7.4162e-05
1.5977	-0.021172	0.001127
0.45627	-0.0010862	-0.0014558
-0.15347	0.016763	0.00066036
0.040141	-0.0062665	-0.00012758
-0.0049951	0.00064387	8.9007e-06];

% Equation (86)
an = cmi(:,1) + T * (cmi(:,2)  + T * cmi(:,3));

% Equation (87)
bn = dmi(:,1) + T * (dmi(:,2) + T * dmi(:,3));

% Equation (85)
x = log(f);

% Equation (83)

alpha = an(1) + x * (an(2) + x * (an(3) + x * (an(4) + x * (an(5) + x * (an(6) + x * (an(7)))))));

% Equation (84)

arg   = bn(1) + x * (bn(2) + x * (bn(3) + x * (bn(4) + x * (bn(5) + x * (bn(6))))));

kappa = exp(arg);

% Equation (82)
gR = kappa * R^(alpha);


return
end