close all; clear; clc;
%%
Data1  	=   load('DataRCS1.dat');
Data2  	=   load('DataRCS2.dat');
%%
figure()
hold on
plot(Data1(:,1),Data1(:,2),'-k','LineWidth',1)
plot(Data2(:,1),Data2(:,2),'--k','LineWidth',1)
hold off
xlabel('$L/\lambda$','Interpret','Latex','FontSize',15)
ylabel('$\sigma_{\theta\theta}/\lambda^2$ [dB]','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
legend('$\theta=30^{\circ}$','$\theta=60^{\circ}$','Interpreter','Latex','FontSize',15,...
       'Location','SouthEast')
xlim([0 4])
ylim([-30 10])
%%
% exportgraphics(gca,'Result.pdf','ContentType','vector')
%%