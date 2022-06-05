function[E_theta,E_phi,E,theta]=FarField_ThetaCut(phi,Data,I)
%%
j       =   sqrt(-1);
Ns      =   1e3;
k       =   2*pi;
%%
phi     =   deg2rad(phi);
theta   =   linspace(-pi,pi,Ns);
%%
[N,~]   =   size(I);
E_phi       =   0;
E_theta     =   0;
for i=1:N-2
    xn      =   Data(i,1);
    yn      =   Data(i,2);
    zn      =   Data(i,3);
    xnm    	=   Data(i,4);
    ynm    	=   Data(i,5);
    znm    	=   Data(i,6);
    xnp    	=   Data(i,7);
    ynp    	=   Data(i,8);
    znp    	=   Data(i,9);
    [lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    k_r     =   k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta));
    dot_p   =   lnp*sin(theta_np)*cos(phi_np)*cos(theta)*cos(phi)+...
                lnp*sin(theta_np)*sin(phi_np)*cos(theta)*sin(phi)-...
                lnp*cos(theta_np)*sin(theta);
    dot_m   =   lnm*sin(theta_nm)*cos(phi_nm)*cos(theta)*cos(phi)+...
                lnm*sin(theta_nm)*sin(phi_nm)*cos(theta)*sin(phi)-...
                lnm*cos(theta_nm)*sin(theta);
    E_theta =   E_theta+I(i+1,1)*(dot_p+dot_m).*exp(j*k_r);
    dot_p   =   -lnp*sin(theta_np)*cos(phi_np)*sin(phi)+...
                lnp*sin(theta_np)*sin(phi_np)*cos(phi);
    dot_m   =   -lnm*sin(theta_nm)*cos(phi_nm)*sin(phi)+...
                lnm*sin(theta_nm)*sin(phi_nm)*cos(phi);
    E_phi  	=   E_phi+I(i+1,1)*(dot_p+dot_m).*exp(j*k_r);
end
%%
E_theta(isnan(E_theta))=0;
E_theta(abs(E_theta)<1e-10)=0;
E_phi(isnan(E_phi))=0;
E_phi(abs(E_phi)<1e-10)=0;
E       =   sqrt(abs(E_theta).^2+abs(E_phi).^2);
end
%%