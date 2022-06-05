function[Data,M]=LShape(dz_,L1,L2)
%%
Ns      =   round(L1/dz_);
if mod(Ns,2)~=0
    Ns      =   Ns+1;
end
dz      =   L1/Ns;
Data1  	=   [];
for i=1:Ns-1
    xn      =   0;
    yn      =   L1-i*dz;
    zn      =   0;
    xnm    	=   0;
    ynm    	=   L1-(i-1)*dz;
    znm    	=   0;
    xnp    	=   0;
    ynp    	=   L1-(i+1)*dz;
    znp    	=   0;
    Data1   =   [Data1;[ xn yn zn xnm ynm znm xnp ynp znp]];
end
%%
Ns      =   round(L2/dz_);
if mod(Ns,2)~=0
    Ns      =   Ns+1;
end
dz      =   L2/Ns;
Data2  	=   [];
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
    Data2   =   [Data2;[ xn yn zn xnm ynm znm xnp ynp znp]];
end
%%
Data    =   [Data1 ; ...
            [Data1(end,7) Data1(end,8) Data1(end,9) ...
            Data1(end,1) Data1(end,2) Data1(end,3) ...
            Data2(1,1) Data2(1,2) Data2(1,3)] ; ...
            Data2];
%%
[Ns,~]  =   size(Data);
M       =   (Ns+1)/2;
end
%%