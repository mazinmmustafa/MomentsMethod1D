function[lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp)
%%
lnm     =   sqrt((xn-xnm)^2+(yn-ynm)^2+(zn-znm)^2);
lnp     =   sqrt((xnp-xn)^2+(ynp-yn)^2+(znp-zn)^2);
%%
theta_nm    =   acos((zn-znm)/lnm);
theta_nm(isnan(theta_nm))=0;
theta_np    =   acos((znp-zn)/lnp);
theta_np(isnan(theta_np))=0;
%%
phi_nm	=   atan((yn-ynm)/(xn-xnm));
phi_nm(isnan(phi_nm))=0;
if (xn-xnm)<0
    phi_nm	=   phi_nm+pi;
end
phi_np	=   atan((ynp-yn)/(xnp-xn));
phi_np(isnan(phi_np))=0;
if (xnp-xn)<0
    phi_np	=   phi_np+pi;
end
end
%%