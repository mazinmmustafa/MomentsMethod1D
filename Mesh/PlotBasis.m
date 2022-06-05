close all; clear; clc;
%%
Data    =   load('Basis.dat');
[N,~]   =   size(Data);
R       =   0.1;
%%
figure('WindowState','maximized')
view([45 45])
% view([0 180])
% xlim([-1 +1]*R)
% ylim([-1 +1]*R)
% zlim([-1 +1]*R)
axis equal
hold on
for i=1:N
    plot3(Data(i,1),Data(i,2),Data(i,3),'.k','MarkerSize',10)
    plot3(Data(i,4),Data(i,5),Data(i,6),'.k','MarkerSize',10)
    plot3(Data(i,7),Data(i,8),Data(i,9),'.k','MarkerSize',10)
end
for i=1:N
%     pause(1E-1);
    plot3([Data(i,1) Data(i,4)],[Data(i,2) Data(i,5)], ...
        [Data(i,3) Data(i,6)],'-r','LineWidth',1)
    plot3([Data(i,4) Data(i,7)],[Data(i,5) Data(i,8)], ...
        [Data(i,6) Data(i,9)],'-b','LineWidth',1)
end
for i=1:N
    if Data(i,10)~=0
        plot3([Data(i,1) Data(i,4)],[Data(i,2) Data(i,5)], ...
        [Data(i,3) Data(i,6)],'-g','LineWidth',1)
        plot3([Data(i,4) Data(i,7)],[Data(i,5) Data(i,8)], ...
        [Data(i,6) Data(i,9)],'-g','LineWidth',1)
        str = sprintf('%d', Data(i,10));
        text(Data(i,4),Data(i,5),Data(i,6),str)
    end
end
hold off
%%






