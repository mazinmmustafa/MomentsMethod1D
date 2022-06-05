close all; clear; clc;
%%
Data    =   load('DataRCS.dat');
%%
figure()
plot(Data(:,1),Data(:,2),'-k','LineWidth',1)
xlabel('$\theta$ [deg]','Interpret','Latex','FontSize',14)
ylabel('$\sigma_{\theta\theta}/\lambda^2$ [dB]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([0 180])
ylim([-30 10])
%%
% exportgraphics(gca,'Result.pdf','ContentType','vector')
%%