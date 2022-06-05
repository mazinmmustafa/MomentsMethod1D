close all; clear; clc;
%%
xn      =   0;
yn      =   0;
zn      =   0;
xnm    	=   0;
ynm    	=   0;
znm    	=   -1;
xnp    	=   0;
ynp    	=   0;
znp    	=   1;
[lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
fprintf('theta_n-\t=\t%0.0f\n',theta_nm*180/pi);
fprintf('phi_n-\t=\t\t%0.0f\n',phi_nm*180/pi);
fprintf('theta_n+\t=\t%0.0f\n',theta_np*180/pi);
fprintf('phi_n+\t=\t\t%0.0f\n',phi_np*180/pi);
figure()
hold on
plot3([xnm xn],[ynm yn],[znm zn],'-b','LineWidth',3)
plot3([xn xnp],[yn ynp],[zn znp],'-r','LineWidth',3)
hold off 
axis equal
view([30 30])
grid on
xlabel('$x$','Interpret','Latex')
ylabel('$y$','Interpret','Latex')
zlabel('$z$','Interpret','Latex')
set(gca,'TickLabel','Latex')
%%
