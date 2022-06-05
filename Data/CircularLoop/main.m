close all; clear; clc;
%%
kOhm    =   1E3;
Data    =   load('Zin.dat');
%%
figure()
hold on
plot(Data(:,1),Data(:,2)/kOhm,'-k','LineWidth',1)
plot(Data(:,1),Data(:,3)/kOhm,'--k','LineWidth',1)
hold off
yticks([-4 -2 0 2 4])
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('$L/\lambda$','Interpret','Latex','FontSize',15)
ylabel('$Z=R+jX$ [$\Omega$]','Interpret','Latex','FontSize',15)
legend('$R$','$X$','Interpreter','Latex','FontSize',15)
ylim([-5 +5])
%%
exportgraphics(gcf,'Figure.pdf','ContentType','vector');
%%