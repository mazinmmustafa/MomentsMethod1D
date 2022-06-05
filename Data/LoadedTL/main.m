close all; clear; clc;
%%
Data    =   load('Zin.dat');
%%
Ns      =   1E3;
c0      =   3E8;
eps0    =   8.8541878128E-12;
mu0     =   4*pi*1E-7;
j       =   sqrt(-1);
%%
L       =   2.5;
S       =   0.03;
a       =   1E-3;
ZL      =   200;
%%
freq    =   logspace(5,8,Ns);
lambda  =   c0./freq;
L_      =   (mu0/pi)*acosh(S/(2*a));
C_      =   (eps0*pi)/acosh(S/(2*a));
Z0      =   sqrt(L_/C_);
Zin     =   Z0*(ZL+j*Z0*tan(2.0*pi*L./lambda))...
            ./(Z0+j*ZL*tan(2.0*pi*L./lambda));
%%
figure()
hold on
plot(Data(:,1),Data(:,2),'-k','LineWidth',1)
plot(Data(:,1),Data(:,3),'--k','LineWidth',1)
plot(freq,real(Zin),'-r','LineWidth',1)
plot(freq,imag(Zin),'--r','LineWidth',1)
plot(Data(:,1),Data(:,2),'-k','LineWidth',1)
plot(Data(:,1),Data(:,3),'--k','LineWidth',1)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('Freq [Hz]','Interpret','Latex','FontSize',15)
ylabel('$Z=R+jX$ [$\Omega$]','Interpret','Latex','FontSize',15)
legend('$R$','$X$',...
    'Interpreter','Latex','FontSize',15,...
    'Location','NorthWest')
set(gca,'XScale','log')
% ylim([-5 +5])
%%
exportgraphics(gcf,'Figure.pdf','ContentType','vector');
%%