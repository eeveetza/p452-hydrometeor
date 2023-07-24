clear all
close all
clc

f = linspace(1,100, 50);
R = 5;
T = [0, 35];

for it = 1:length(T)

    for i = 1:length(f)


        gammaR(i) =  kappa_rain(1, 0, f(i), T(it), R);

    end

    loglog(f, gammaR)
    set(gca,'YLim', [1e-4, 1e1]);
    grid on
    hold on

    set(gca,'FontSize',14)

end
xlabel('Frequency (GHz)')
ylabel('Specific attenuation (dB/km)')
legend('T = 0 deg C', 'T = 35 deg C')


figure

R = 100;
T = [0, 35];

for it = 1:length(T)

    for i = 1:length(f)


        gammaR(i) =  kappa_rain(1, 0, f(i), T(it), R);

    end

    loglog(f, gammaR)
    set(gca,'YLim', [1e-3, 1e2]);
    grid on
    hold on

    set(gca,'FontSize',14)

end
xlabel('Frequency (GHz)')
ylabel('Specific attenuation (dB/km)')
legend('T = 0 deg C', 'T = 35 deg C')