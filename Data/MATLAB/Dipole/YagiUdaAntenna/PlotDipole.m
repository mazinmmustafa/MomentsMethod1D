function[]=PlotDipole(N,M,Data)
%% 
figure()
hold on
for n=1:N
    xn      =   Data(n,1);   
    zn      =   Data(n,3);
    xnm     =   Data(n,4);   
    znm  	=   Data(n,6);
    xnp     =   Data(n,7);   
    znp  	=   Data(n,9);
    plot([xnm xn],[znm zn],'-','LineWidth',1.5)
    plot([xn xnp],[zn znp],'-','LineWidth',1.5)
    if n==M
        plot([xnm xn],[znm zn],'-r','LineWidth',3)
        plot([xn xnp],[zn znp],'-r','LineWidth',3)
    end
    axis equal
    xlim([-1 +1]/4)
    ylim([-1 +1]/4)
    grid on
    xlabel('$x/\lambda$','Interpret','Latex','FontSize',14)
    ylabel('$z/\lambda$','Interpret','Latex','FontSize',14)
    set(gca,'TickLabel','Latex','FontSize',14)
end
hold off
end
%%