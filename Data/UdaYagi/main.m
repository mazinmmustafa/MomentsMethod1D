close all; clear; clc;
%%
Data    =   load('In.dat');
N       =   length(Data);
Element =   1:N;
%%
figure()
hold on
plot(Element,Data(:,1)/max(Data(:,1)),'-k','LineWidth',1)
plot(Element,Data(:,1)/max(Data(:,1)),'.k','MarkerSize',12)
hold off
set(gca,'TickLabel','Latex','FontSize',15)
xlabel('Element','Interpret','Latex','FontSize',15)
ylabel('$|I_{n}|$','Interpret','Latex','FontSize',15)
ylim([0 1])
%%
exportgraphics(gcf,'FigureCurrent.pdf','ContentType','vector');
%%