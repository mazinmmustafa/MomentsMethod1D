function[]=PlotCurrent(N,Data,I)
%%
mA      =   1e-3;
I_      =   [0;I;0];
figure()
hold on
for n=1:N
    zn      =   Data(n,3); 
    znm  	=   Data(n,6);
    plot([znm zn],real([I_(n) I_(n+1)])/mA,'-k','LineWidth',1)
    plot([znm zn],imag([I_(n) I_(n+1)])/mA,'--k','LineWidth',1)
end
znp  	=   Data(n,9);
plot([zn znp],real([I_(n+1) I_(n+2)])/mA,'-k','LineWidth',1)
plot([zn znp],imag([I_(n+1) I_(n+2)])/mA,'--k','LineWidth',1)
hold off
xlabel('$l/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$I_{n}(l)$ [mA]','Interpret','Latex','FontSize',14)
legend('$\mathcal{R}e$','$\mathcal{I}m$','Interpreter','Latex','Location',...
        'NorthEast','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
end
%%