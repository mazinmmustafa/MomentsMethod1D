close all; clear; clc;
%%
Data    =   load('DataPeakCurrent.dat');
%%
figure()
plot(0:length(Data)-1,Data(:,1),'-k','LineWidth',1)
hold on
plot(0:length(Data)-1,Data(:,1),'.k','MarkerSize',14)
hold off
xlabel('Element','Interpret','Latex','FontSize',15)
ylabel('$|I_{n}|$','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
xlim([0 length(Data)-1])
ylim([0 1])
%%
exportgraphics(gca,'Result.pdf','ContentType','vector')
%%