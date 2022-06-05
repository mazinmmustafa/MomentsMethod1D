function[Data]=VerticalDipoleS(dz,L)
%%
Ns      =   round(L/dz,0);
%%
dz      =   L/Ns;
Data    =   zeros(Ns-1,9);
%%
for i=1:Ns-1
    xn      =   0;
    yn      =   0;
    zn      =   i*dz;
    xnm    	=   0;
    ynm    	=   0;
    znm    	=   (i-1)*dz;
    xnp    	=   0;
    ynp    	=   0;
    znp    	=   (i+1)*dz;
    Data(i,:)   =   [xn yn zn xnm ynm znm xnp ynp znp];
end
%%
Data(:,3)   =   Data(:,3)-L/2;  
Data(:,6)   =   Data(:,6)-L/2;
Data(:,9)   =   Data(:,9)-L/2;
end
%%