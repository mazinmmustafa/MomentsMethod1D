close all; clear; clc;
%%
Ns      =   201;
L_min   =   0.001;
L_max   =   2;
%%
L       =   linspace(L_min,L_max,Ns);
a       =   L/(2*74.2);
%%
N_      =   zeros(1,Ns);
Yin     =   zeros(1,Ns);
count	=   1;
tic;
for i=1:Ns
    N       =   round(15*L(i));
    if mod(N,2)==0
        N       =   N+1;
    end
    N_(1,i) =   N;
	fprintf("Step:\t%3d/100,\tN\t=\t%2d,\tL\t=\t%0.3f\n",round(100*count/Ns),N,L(i));
    count   =   count+1;
    [Data,N,M]=Dipole(N,L(i));
    I=MM(N,M,Data,a(i));
    Yin(1,i)=   I(M);
    clc;
end
toc;
%%
mS      =   1e-3;
figure()
hold on
plot(L,real(Yin)/mS,'-k','LineWidth',1)
plot(L,imag(Yin)/mS,'--k','LineWidth',1)
hold off
xlabel('$L/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$Y_{\textrm{in}}$ [mS]','Interpret','Latex','FontSize',14)
legend('$G$','$B$','Interpreter','Latex','Location',...
        'NorthEast','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([0 L_max])
%%
figure()
stem(L,N_,'-k')
xlabel('$L/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$N$','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([0 L_max])
%%
% exportgraphics(gca,'Admittance.pdf','ContentType','vector');
%%
