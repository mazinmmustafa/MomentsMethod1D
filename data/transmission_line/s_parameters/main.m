close all; clear; clc;

Ns = 1000;
MHz = 1E6;
cm = 1E-2;
mm = 1E-3;
eta0 = 120*pi;
c0 = 3E8;

L = 250*cm;
h = 3*cm;
a = 1*mm;
Z0 = 50;

freq_min = 0.1*MHz;
freq_max = 100*MHz;
freq = linspace(freq_min, freq_max, Ns);
lambda = c0./freq;
Zc = (1.0/pi)*eta0*acosh(h/(2.0*a));

j = sqrt(-1);
beta_l = 2*pi*L./lambda;
Gamma_s = (Zc-Z0)./(Zc+Z0);
Gamma_l = (Z0-Zc)./(Z0+Zc);
Gamma_0 = Gamma_l*exp(-j*2*beta_l);
S11 = (Gamma_s+Gamma_0)./(1+Gamma_s*Gamma_0);
S21 = (1+Gamma_s).*(1+Gamma_l).*exp(j*beta_l)./(1+Gamma_s*Gamma_0);

figure()
hold on
plot(freq/MHz,20*log10(abs(S11)))
plot(freq/MHz,20*log10(abs(S21)))
hold off
ylim([-30 0])