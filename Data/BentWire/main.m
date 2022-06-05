close all; clear; clc;
%%
GHz     =   1E9;
j       =   sqrt(-1);
Data    =   load('Zin.dat');
Freq    =   Data(:,1);
Zin     =   Data(:,2)+j*Data(:,3);
Z0      =   50;
S11     =   (Zin-Z0)./(Zin+Z0);
%%
DataHFSS 	=   csvread('HFSS.csv');
% DataHFSS 	=   csvread('HFSS_S.csv');
%%
figure()
hold on
plot(Freq/GHz,20*log10(abs(S11)),'-k','LineWidth',1)
plot(DataHFSS(:,1),DataHFSS(:,2),'--k','LineWidth',1)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('Freq [GHz]','Interpret','Latex','FontSize',15)
ylabel('$|S_{11}|$ [dB]','Interpret','Latex','FontSize',15)
legend('MM','HFSS','Interpreter','Latex','FontSize',15,...
        'Location','SouthEast')
xlim([0 2.0])
ylim([-20 0])
% pbaspect([2 1 1])
%%
exportgraphics(gcf,'Figure.pdf','ContentType','vector');
%%