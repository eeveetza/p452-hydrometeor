f = 20; 
R = 100;
T = [0, 35];

phi = linspace(0,pi,1000);

for it = 1:length(T)

for i = 1:length(phi)


eta(i) = raindrop_bscross_section(phi(i), f, R, T(it));

end 

plot(phi*180/pi, 10*log10(eta))
set(gca,'YLim', [-32, -28]);
grid on
hold on

set(gca,'FontSize',14)

end
xlabel('Scattering angle (Degrees)')
ylabel('Bi-static cross section (dB)')
legend('T = 0 deg C', 'T = 35 deg C')


figure


f = 100; 
R = 100;
T = [0, 35];

phi = linspace(0,pi,1000);

for it = 1:length(T)

for i = 1:length(phi)


eta(i) = raindrop_bscross_section(phi(i), f, R, T(it));

end 

plot(phi*180/pi, 10*log10(eta))
set(gca,'YLim', [-30, -14]);
grid on
hold on

set(gca,'FontSize',14)

end
xlabel('Scattering angle (Degrees)')
ylabel('Bi-static cross section (dB)')
legend('T = 0 deg C', 'T = 35 deg C')