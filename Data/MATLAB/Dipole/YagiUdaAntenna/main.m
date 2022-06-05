close all; clear; clc;
%%
L       =   2;
a       =   L/200;  
theta_i =   45;
phi_i   =   0;
E_TM    =   1;
E_TE    =   0;
N       =   31;
%%
theta_i =   theta_i*pi/180;
phi_i   =   phi_i*pi/180;
%%
tic;
[Data,N]=Dipole(N,L);
I=MM(N,Data,a,E_TM,E_TE,theta_i,phi_i);
toc;
% PlotCurrent(N,Data,I);
% plot(linspace(-L/2,L/2,N),abs(I)/1e-3)
%%
Ns      =   1e2;
theta   =   linspace(0,pi,Ns);
phi     =   0;
eta     =   120*pi;
sum     =   0;
for n=1:N
    xn      =   Data(n,1);   
    yn      =   Data(n,2);
    zn      =   Data(n,3);
    xnm     =   Data(n,4);   
    ynm  	=   Data(n,5);
    znm  	=   Data(n,6);
    xnp     =   Data(n,7);   
    ynp     =   Data(n,8);   
    znp  	=   Data(n,9);
    dot_p   =   (xnp-xn)*cos(theta)*cos(phi)+(ynp-yn)*cos(theta)*sin(phi)-(znp-zn)*sin(theta);
    dot_m   =   (xn-xnm)*cos(theta)*cos(phi)+(yn-ynm)*cos(theta)*sin(phi)-(zn-znm)*sin(theta);
    sum     =   sum+(dot_p.*kn_p(Data,theta,phi,n)+dot_m.*kn_m(Data,theta,phi,n))*I(n);
end
sigma   =   4*pi*abs(0.5*eta*sum).^2;
%%
figure()
plot(theta*180/pi,10*log10(sigma))
ylim([-20 10])
%%
function[kn]=kn_p(Data,theta,phi,n)
j       =   sqrt(-1);
k       =   2*pi;
xn      =   Data(n,1);   
yn      =   Data(n,2);
zn      =   Data(n,3);
xnp     =   Data(n,7);   
ynp     =   Data(n,8);   
znp  	=   Data(n,9);
%%
L_r     =   (xnp-xn)*sin(theta)*cos(phi)+(ynp-yn)*sin(theta)*sin(phi)+(znp-zn)*cos(theta);
if L_r==0
    kn    	=   exp(j*k*(xnp*sin(theta)*cos(phi)+ynp*sin(theta)*sin(phi)+znp*cos(theta)))/2;
else
    factor  =   exp(j*k*(xnp*sin(theta)*cos(phi)+ynp*sin(theta)*sin(phi)+znp*cos(theta)));
    kn    	=   factor.*(exp(-j*k*L_r).*(1+j*k*L_r)-1)./(k^2*L_r.^2);
end
end
%%
function[kn]=kn_m(Data,theta,phi,n)
j       =   sqrt(-1);
k       =   2*pi;
xn      =   Data(n,1);   
yn      =   Data(n,2);
zn      =   Data(n,3);
xnm     =   Data(n,4);   
ynm  	=   Data(n,5);
znm  	=   Data(n,6);
%%
L_r     =   (xn-xnm)*sin(theta)*cos(phi)+(yn-ynm)*sin(theta)*sin(phi)+(zn-znm)*cos(theta);
if L_r==0
    kn    	=   exp(j*k*(xnm*sin(theta)*cos(phi)+ynm*sin(theta)*sin(phi)+znm*cos(theta)))/2;
else
    factor  =   exp(j*k*(xnm*sin(theta)*cos(phi)+ynm*sin(theta)*sin(phi)+znm*cos(theta)));
    kn    	=   factor.*(exp(+j*k*L_r).*(1-j*k*L_r)-1)./(k^2*L_r.^2);
end
end
%%
