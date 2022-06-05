close all; clear; clc;
%%
mA      =   1E-3;
Data    =   load('In.dat');
L       =   linspace(0,1.5,length(Data));
%%
figure()
hold on
plot(L,Data(:,1)/mA,'-k','LineWidth',1)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('$L/\lambda$','Interpret','Latex','FontSize',15)
ylabel('$|I_{n}|$ [mA]','Interpret','Latex','FontSize',15)
ylim([0 10])
%%
exportgraphics(gcf,'Figure.pdf','ContentType','vector');
%%