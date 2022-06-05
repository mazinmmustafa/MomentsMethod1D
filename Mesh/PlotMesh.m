close all; clear; clc;
%%
Data    =   load('Mesh.dat');
[N,~]   =   size(Data);
R       =   2;
%%
figure()
view([45 45])
axis equal
hold on
for i=1:N
    plot3(Data(i,1),Data(i,2),Data(i,3),'.k','MarkerSize',10)
    plot3(Data(i,4),Data(i,5),Data(i,6),'.k','MarkerSize',10)
    if Data(i,7)~=0
        plot3([Data(i,1) Data(i,4)],[Data(i,2) Data(i,5)], ...
        [Data(i,3) Data(i,6)],'-r','LineWidth',1)
    else
        plot3([Data(i,1) Data(i,4)],[Data(i,2) Data(i,5)], ...
        [Data(i,3) Data(i,6)],'-b','LineWidth',1)
    end
%     str = sprintf('%d', i);
%     text(0.5*(Data(i,1)+Data(i,4)),0.5*(Data(i,2)+Data(i,5)),...
%         0.5*(Data(i,3)+Data(i,6)),str)
end
hold off
% xlim([-1 +1]*R)
% ylim([-1 +1]*R)
% zlim([-1 +1]*R)
%%






