close all; clear; clc;
%%
c0      =   3e8;
MHz     =   1e6;
%%
f       =   900*MHz;
lambda  =   c0/f;
L       =   1/lambda;
L1     	=   0.7*L;
L2     	=   0.3*L;
a       =   (L1+L2)/200;
E_TE    =   0;
E_TM    =   1;
theta_i =   30;
phi_i   =   30;
%%
dz      =   1/15;
%%
Data    =   LShape(dz,L1,L2);
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
    plot3([xnm xn],[ynm yn],[znm zn],'-','LineWidth',1)
    plot3([xn xnp],[yn ynp],[zn znp],'-','LineWidth',1)
    hold off
    
end
view([30 30])
axis equal
grid on
xlabel('$x$','Interpret','Latex')
ylabel('$y$','Interpret','Latex')
zlabel('$z$','Interpret','Latex')
set(gca,'TickLabel','Latex')
%%
tic;
I       =   MM_PlaneWave(Data,Nss,theta_i,phi_i,E_TE,E_TM,a);
I       =   [0;I;0];
%%
L_    	=   zeros(1,Nss);
xn      =   Data(1,1);
yn      =   Data(1,2);
zn      =   Data(1,3);
xnm    	=   Data(1,4);
ynm    	=   Data(1,5);
znm    	=   Data(1,6);
xnp    	=   Data(1,7);
ynp    	=   Data(1,8);
znp    	=   Data(1,9);
[lnm,~,~,~,~,~]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
sum     =   lnm;
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
    [~,lnp,~,~,~,~]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    L_(1,i) =   sum;
    sum     =   sum+lnp;
end
L_      =   [0 L_ sum];
%%
V           =   zeros(Nss,1);
for m=1:Nss
    xm      =   Data(m,1);
    ym      =   Data(m,2);
    zm      =   Data(m,3);
    xmm    	=   Data(m,4);
    ymm    	=   Data(m,5);
    zmm    	=   Data(m,6);
    xmp    	=   Data(m,7);
    ymp    	=   Data(m,8);
    zmp    	=   Data(m,9);
    [lmm,lmp,theta_mm,phi_mm,theta_mp,phi_mp]=rn(xm,ym,zm,xmm,ymm,zmm,xmp,ymp,zmp);
 	l_mx        =   sin((theta_mp+theta_mm)/2)*cos((phi_mp+phi_mm)/2);
  	l_my        =   sin((theta_mp+theta_mm)/2)*sin((phi_mp+phi_mm)/2);
   	l_mz        =   cos((theta_mp+theta_mm)/2);
    [Ex,Ey,Ez]=Ei(xm,ym,zm,theta_i,phi_i,E_TM,E_TE);
    V(m,1)      =   l_mx*Ex+l_my*Ey+l_mz*Ez;
end
V           =   max(abs(V));
%%
mA      =   1e-3;
figure()
hold on
plot(L_,real(I)/mA,'-k','LineWidth',1)
plot(L_,imag(I)/mA,'--k','LineWidth',1)
hold off
xlabel('$l/\lambda$','Interpret','Latex')
ylabel('$I_{n}(l)$ [mA]','Interpret','Latex')
legend('Re','Im','Interpreter','Latex','Location','NorthEast')
set(gca,'TickLabel','Latex')
%%
% exportgraphics(gca,'Currents.pdf','ContentType','vector');
%%
mA      =   1e-3;
[N,~]   =   size(I);
figure()
plot(L_,abs(I)/mA,'-ok','MarkerFace','k','MarkerSize',2,'LineWidth',1)
% plot(L_-L/2,abs(I)/V,'-ok','MarkerFace','k','MarkerSize',2,'LineWidth',1)
hold on
for i=2:N-1
    plot([L_(i-1) L_(i)],[0 abs(I(i))]/mA,'-k','LineWidth',0.1)
    plot([L_(i) L_(i+1)],[abs(I(i)) 0]/mA,'-k','LineWidth',0.1)
%     plot([L_(i-1) L_(i)],[0 abs(I(i))]/V,'-k','LineWidth',0.1)
%     plot([L_(i) L_(i+1)],[abs(I(i)) 0]/V,'-k','LineWidth',0.1)
end
hold off
xlabel('$l/\lambda$','Interpret','Latex')
ylabel('$|I_{n}(l)|$ [mA]','Interpret','Latex')
set(gca,'TickLabel','Latex')
%%
% exportgraphics(gca,'CurrentsAbs.pdf','ContentType','vector');
%%
[E_theta,~,~,phi]=FarField_PhiCut(90,Data,I);
E_theta 	=   E_theta./max(abs(E_theta));
figure()
PolardB(phi,20*log10(abs(E_theta)),[-25 0],6,'-k')
title('$20\log_{10}|E_{\theta}|$ [dB]','Interpret','Latex')
pax = gca;
pax.RAxisLocation   =   90;
pax.ThetaDir = 'CounterClockwise';
pax.ThetaZeroLocation = 'right';
%%
% exportgraphics(gca,'FarFieldAzimuth.pdf','ContentType','vector');
%%
[E_theta,~,~,theta]=FarField_ThetaCut(0,Data,I);
E_theta    	=   E_theta./max(abs(E_theta));
figure()
PolardB(theta,20*log10(abs(E_theta)),[-25 0],6,'-k')
title('$20\log_{10}|E_{\theta}|$ [dB]','Interpret','Latex')
%%
% exportgraphics(gca,'FarFieldElevation.pdf','ContentType','vector');
%%
[sigma_theta,sigma_phi,theta]=Sigma_theta(Data,I,60);
%%
figure()
plot(theta*180/pi,10*log10(abs(sigma_theta)),'-k','LineWidth',1)
xlabel('$\theta$ [deg]','Interpret','Latex')
ylabel('$\sigma_{\theta\theta}/\lambda^2 $ [dB]','Interpret','Latex')
set(gca,'TickLabel','Latex')
xlim([-180 180])
% ylim([-30 10])
%%
% exportgraphics(gca,'RCS_theta_1.pdf','ContentType','vector');
%%
figure()
plot(theta*180/pi,10*log10(abs(sigma_phi)),'-k','LineWidth',1)
xlabel('$\theta$ [deg]','Interpret','Latex')
ylabel('$\sigma_{\varphi\theta}/\lambda^2 $ [dB]','Interpret','Latex')
set(gca,'TickLabel','Latex')
xlim([-180 180])
% ylim([-30 10])
%%
% exportgraphics(gca,'RCS_phi_1.pdf','ContentType','vector');
%%
[sigma_theta,sigma_phi,phi]=Sigma_phi(Data,I,60);
%%
figure()
plot(phi*180/pi,10*log10(abs(sigma_theta)),'-k','LineWidth',1)
xlabel('$\varphi$ [deg]','Interpret','Latex')
ylabel('$\sigma_{\theta\theta}/\lambda^2 $ [dB]','Interpret','Latex')
set(gca,'TickLabel','Latex')
xlim([-180 180])
% ylim([-30 10])
%%
% exportgraphics(gca,'RCS_theta_2.pdf','ContentType','vector');
%%
figure()
plot(phi*180/pi,10*log10(abs(sigma_phi)),'-k','LineWidth',1)
xlabel('$\varphi$ [deg]','Interpret','Latex')
ylabel('$\sigma_{\varphi\theta}/\lambda^2 $ [dB]','Interpret','Latex')
set(gca,'TickLabel','Latex','FontSize',15)
xlim([-90 90])
ylim([-30 10])
%%
% exportgraphics(gca,'RCS_phi_2.pdf','ContentType','vector');
%% 
function[Ex,Ey,Ez]=Ei(x,y,z,theta_i,phi_i,E_TM,E_TE)
j           =   sqrt(-1);
k           =   2*pi;
theta_i     =   deg2rad(theta_i);
phi_i       =   deg2rad(phi_i);
k_r         =   k*(x*sin(theta_i)*cos(phi_i)+y*sin(theta_i)*sin(phi_i)+z*cos(theta_i));
Ex          =   (E_TM*(cos(theta_i)*cos(phi_i))-E_TE*(sin(phi_i)))*exp(j*k_r);
Ey          =   (E_TM*(cos(theta_i)*sin(phi_i))+E_TE*(cos(phi_i)))*exp(j*k_r);
Ez          =   (-E_TM*sin(theta_i))*exp(j*k_r);
end
%%
function[sigma_theta,sigma_phi,theta]=Sigma_theta(Data,I,phi)
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
theta   =   linspace(-pi,pi,1e3);
phi   	=   phi*pi/180;
sigma_theta     =   0;
sigma_phi       =   0;
[N,~]   = size(I);
for n=1:N-2
    xn      =   Data(n,1);
    yn      =   Data(n,2);
    zn      =   Data(n,3);
    xnm    	=   Data(n,4);
    ynm    	=   Data(n,5);
    znm    	=   Data(n,6);
    xnp    	=   Data(n,7);
    ynp    	=   Data(n,8);
    znp    	=   Data(n,9);
    [lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    dot1_np =   sin(theta_np)*cos(phi_np)*cos(theta)*cos(phi)+...
                sin(theta_np)*sin(phi_np)*cos(theta)*sin(phi)-...
                cos(theta_np)*sin(theta);
    dot1_nm =   sin(theta_nm)*cos(phi_nm)*cos(theta)*cos(phi)+...
                sin(theta_nm)*sin(phi_nm)*cos(theta)*sin(phi)-...
                cos(theta_nm)*sin(theta);
    dot2_np =   -sin(theta_np)*cos(phi_np)*sin(phi)+...
                sin(theta_np)*sin(phi_np)*cos(phi);
    dot2_nm =   -sin(theta_nm)*cos(phi_nm)*sin(phi)+...
                sin(theta_nm)*sin(phi_nm)*cos(phi);
    Term1   =   (eta/4)*I(n+1,1)*(lnp*dot1_np+lnm*dot1_nm).*...
        exp(j*k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta)));
    Term2   =   (eta/4)*I(n+1,1)*(lnp*dot2_np+lnm*dot2_nm).*...
        exp(j*k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta)));
    sigma_theta     =   sigma_theta+Term1;
    sigma_phi       =   sigma_phi+Term2;
end
sigma_theta	=   4*pi*abs(sigma_theta).^2;
sigma_phi 	=   4*pi*abs(sigma_phi).^2;
end
%%
function[sigma_theta,sigma_phi,phi]=Sigma_phi(Data,I,theta)
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
theta   =   theta*pi/180;
phi   	=   linspace(-pi,pi,1e3);
sigma_theta     =   0;
sigma_phi       =   0;
[N,~]   = size(I);
for n=1:N-2
    xn      =   Data(n,1);
    yn      =   Data(n,2);
    zn      =   Data(n,3);
    xnm    	=   Data(n,4);
    ynm    	=   Data(n,5);
    znm    	=   Data(n,6);
    xnp    	=   Data(n,7);
    ynp    	=   Data(n,8);
    znp    	=   Data(n,9);
    [lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    dot1_np =   sin(theta_np)*cos(phi_np)*cos(theta)*cos(phi)+...
                sin(theta_np)*sin(phi_np)*cos(theta)*sin(phi)-...
                cos(theta_np)*sin(theta);
    dot1_nm =   sin(theta_nm)*cos(phi_nm)*cos(theta)*cos(phi)+...
                sin(theta_nm)*sin(phi_nm)*cos(theta)*sin(phi)-...
                cos(theta_nm)*sin(theta);
    dot2_np =   -sin(theta_np)*cos(phi_np)*sin(phi)+...
                sin(theta_np)*sin(phi_np)*cos(phi);
    dot2_nm =   -sin(theta_nm)*cos(phi_nm)*sin(phi)+...
                sin(theta_nm)*sin(phi_nm)*cos(phi);
    Term1   =   (eta/4)*I(n+1,1)*(lnp*dot1_np+lnm*dot1_nm).*...
        exp(j*k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta)));
    Term2   =   (eta/4)*I(n+1,1)*(lnp*dot2_np+lnm*dot2_nm).*...
        exp(j*k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta)));
    sigma_theta     =   sigma_theta+Term1;
    sigma_phi       =   sigma_phi+Term2;
end
sigma_theta	=   4*pi*abs(sigma_theta).^2;
sigma_phi 	=   4*pi*abs(sigma_phi).^2;
end
%%

