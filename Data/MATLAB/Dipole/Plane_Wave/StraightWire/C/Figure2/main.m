close all; clear; clc;
%%
Data  	=   load('DataRCS.dat');
%%
figure()
plot(Data(:,1),Data(:,2),'-k','LineWidth',1)
xlabel('$L/\lambda$','Interpret','Latex','FontSize',15)
ylabel('$\sigma_{\theta\theta}/\lambda^2$ [dB]','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
xlim([0 4])
ylim([-30 10])
%%
% exportgraphics(gca,'Result.pdf','ContentType','vector')
%%