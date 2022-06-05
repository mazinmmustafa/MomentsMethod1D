close all; clear; clc;
%%
mS      =   1E-3;
Data    =   load('Yin.dat');
%%
figure()
hold on
plot(Data(:,1),Data(:,2)/mS,'-k','LineWidth',1)
plot(Data(:,1),Data(:,3)/mS,'--k','LineWidth',1)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('$L/\lambda$','Interpret','Latex','FontSize',15)
ylabel('$Y=G+jB$ [mS]','Interpret','Latex','FontSize',15)
legend('$G$','$B$','Interpreter','Latex','FontSize',15)
%%
exportgraphics(gcf,'Figure.pdf','ContentType','vector');
%%