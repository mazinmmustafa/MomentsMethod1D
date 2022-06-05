function[Data,N]=Dipole(N,L)
%%
if mod(N,2)==0
    N       =   N+1;
end
Data    =   zeros(N,9);
dl      =   L/(N+1);
for n=1:N
    xn      =   0;   
    yn      =   0;
    zn      =   n*dl-L/2;
    xnm     =   0;   
    ynm 	=   0;
    znm  	=   (n-1)*dl-L/2;
    xnp     =   0;   
    ynp 	=   0;
    znp  	=   (n+1)*dl-L/2;
    Data(n,:)   =   [ xn yn zn xnm ynm znm xnp ynp znp ];
end
end
%%