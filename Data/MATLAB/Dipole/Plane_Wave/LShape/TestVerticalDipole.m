close all; clear; clc;
%%
Ns      =   10;
%%
L       =   0.5;
a       =   0.005;
%%
[Data,M]=VerticalDipole(Ns,L);
[Nss,~] =   size(Data);
%%
figure()
for i=1:Nss
    xn      =   Data(i,1);
    yn      =   Data(i,2);
    zn      =   Data(i,3);
    xnm    	=   Data(i,4);
    ynm    	=   Data(i,5);
    znm    	=   Data(i,6);
    xnp    	=   Data(i,7);
    ynp    	=   Data(i,8);
    znp    	=   Data(i,9);
    hold on
    if i==M
        plot3([xnm xn],[ynm yn],[znm zn],'-r','LineWidth',3)
        plot3([xn xnp],[yn ynp],[zn znp],'-r','LineWidth',3)
    else
        plot3([xnm xn],[ynm yn],[znm zn],'-','LineWidth',1)
        plot3([xn xnp],[yn ynp],[zn znp],'-','LineWidth',1)
    end
    hold off
end
axis equal
view([30 30])
grid on
xlabel('$x$','Interpret','Latex')
ylabel('$y$','Interpret','Latex')
zlabel('$z$','Interpret','Latex')
set(gca,'TickLabel','Latex')
%%