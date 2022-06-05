close all; clear; clc;
%%
MHz     =   1E6;
j       =   sqrt(-1);
Data    =   load('Zin.dat');
Freq    =   Data(:,1);
Zin     =   Data(:,2)+j*Data(:,3);
Z0      =   50;
S11     =   (Zin-Z0)./(Zin+Z0);
%%
figure()
hold on
plot(Freq/MHz,20*log10(abs(S11)),'-xk','LineWidth',1)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('Freq [MHz]','Interpret','Latex','FontSize',15)
ylabel('$|S_{11}|$ [dB]','Interpret','Latex','FontSize',15)
xlim([25 300])
ylim([-20 0])
% pbaspect([2 1 1])
%%
exportgraphics(gcf,'Figure3.pdf','ContentType','vector');
%%
figure()
hold on
plot(Freq/MHz,real(Zin),'-xk','LineWidth',1)
plot(Freq/MHz,imag(Zin),'--xk','LineWidth',1)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('Freq [MHz]','Interpret','Latex','FontSize',15)
ylabel('$Z=R+jX$ [$\Omega$]','Interpret','Latex','FontSize',15)
legend('$R$','$X$','Interpreter','Latex','FontSize',15,...
    'Location','SouthEast')
xlim([25 300])
ylim([-300 400])
% pbaspect([2 1 1])
%%
exportgraphics(gcf,'Figure4.pdf','ContentType','vector');
%%