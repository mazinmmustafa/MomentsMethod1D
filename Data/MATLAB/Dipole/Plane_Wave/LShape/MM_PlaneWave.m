function[In]=MM_PlaneWave(Data,Nss,theta_i,phi_i,E_TE,E_TM,a)
%%
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
theta_i =   theta_i*pi/180;
phi_i   =   phi_i*pi/180;
%%
Zmn   	=   zeros(Nss,Nss);
Vm   	=   zeros(Nss,1);
for m=1:Nss
    for n=1:Nss
        %% Data m
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
        %% Data n
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
        %%
        if m==n
            Zmn(m,n)    =   SingularTerm(lmp)+SingularTerm(lmm)+Term_pm+Term_mp;
        elseif m==n+1 
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+SingularTerm2(lmm); 
        elseif n==m+1 
            Zmn(m,n)    =   Term_pp+Term_mm+Term_mp+SingularTerm2(lmp);
        else
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_mp;  
        end  
    end
    ki_l_mp =   k*(sin(theta_i)*cos(phi_i)*sin(theta_mp)*cos(phi_mp)+...
                sin(theta_i)*sin(phi_i)*sin(theta_mp)*sin(phi_mp)+...
                cos(theta_i)*cos(phi_mp));
    ki_l_mm =   k*(sin(theta_i)*cos(phi_i)*sin(theta_mm)*cos(phi_mm)+...
                sin(theta_i)*sin(phi_i)*sin(theta_mm)*sin(phi_mm)+...
                cos(theta_i)*cos(phi_mm));
	dot1_mp =   sin(theta_mp)*cos(phi_mp)*cos(theta_i)*cos(phi_i)+...
                sin(theta_mp)*sin(phi_mp)*cos(theta_i)*sin(phi_i)-...
                cos(theta_mp)*sin(theta_i);
    dot1_mm =   sin(theta_mm)*cos(phi_mm)*cos(theta_i)*cos(phi_i)+...
                sin(theta_mm)*sin(phi_mm)*cos(theta_i)*sin(phi_i)-...
                cos(theta_mm)*sin(theta_i);
    dot2_mp =   -sin(theta_mp)*cos(phi_mp)*sin(phi_i)+...
                sin(theta_mp)*sin(phi_mp)*cos(phi_i);
    dot2_mm =   -sin(theta_mm)*cos(phi_mm)*sin(phi_i)+...
                sin(theta_mm)*sin(phi_mm)*cos(phi_i);
    Xmp     =   @(l) (1-l/lmp).*exp(+j*ki_l_mp*l);
  	Xmm     =   @(l) (1-l/lmm).*exp(-j*ki_l_mm*l);
   	Vm(m,1)	=   exp(j*k*(xm*sin(theta_i)*cos(phi_i)+ym*sin(theta_i)*sin(phi_i)+zm*cos(theta_i)))*...
                (E_TM*dot1_mp*Quad(Xmp,0,lmp)+E_TM*dot1_mm*Quad(Xmm,0,lmm)+...
                E_TE*dot2_mp*Quad(Xmp,0,lmp)+E_TE*dot2_mm*Quad(Xmm,0,lmm));
end
%%
In      =   Zmn\Vm;
%% Term ++
function[Termpp]=Term_pp()
dot_pp  =   sin(theta_mp)*cos(phi_mp)*sin(theta_np)*cos(phi_np)+...
            sin(theta_mp)*sin(phi_mp)*sin(theta_np)*sin(phi_np)+...
          	cos(theta_mp)*cos(theta_np);
R_pp    =   @(l,l_) Rmn_pp(l,l_,xm,ym,zm,xn,yn,zn,theta_mp,phi_mp,theta_np,phi_np);
g_pp  	=   @(l,l_) exp(-j*k*R_pp(l,l_))./(4*pi*R_pp(l,l_));
Z1_pp   =   @(l,l_) j*k*eta*dot_pp*g_pp(l,l_).*(1-l/lmp).*(1-l_/lnp);
Z2_pp   =   @(l,l_) -j*(eta/k)*(1/(lmp*lnp))*g_pp(l,l_);
Termpp 	=   Quad2(Z1_pp,0,lmp,0,lnp)+Quad2(Z2_pp,0,lmp,0,lnp);
end
%% Term --
function[Termmm]=Term_mm()
dot_mm  =   sin(theta_mm)*cos(phi_mm)*sin(theta_nm)*cos(phi_nm)+...
            sin(theta_mm)*sin(phi_mm)*sin(theta_nm)*sin(phi_nm)+...
            cos(theta_mm)*cos(theta_nm);
R_mm    =   @(l,l_) Rmn_mm(l,l_,xm,ym,zm,xn,yn,zn,theta_mm,phi_mm,theta_nm,phi_nm);
g_mm  	=   @(l,l_) exp(-j*k*R_mm(l,l_))./(4*pi*R_mm(l,l_));
Z1_mm   =   @(l,l_) j*k*eta*dot_mm*g_mm(l,l_).*(1-l/lmm).*(1-l_/lnm);
Z2_mm   =   @(l,l_) -j*(eta/k)*(1/(lmm*lnm))*g_mm(l,l_);
Termmm 	=   Quad2(Z1_mm,0,lmm,0,lnm)+Quad2(Z2_mm,0,lmm,0,lnm);
end
%% Term +-
function[Termpm]=Term_pm()
dot_pm  =   sin(theta_mp)*cos(phi_mp)*sin(theta_nm)*cos(phi_nm)+...
            sin(theta_mp)*sin(phi_mp)*sin(theta_nm)*sin(phi_nm)+...
            cos(theta_mp)*cos(theta_nm);
R_pm    =   @(l,l_) Rmn_pm(l,l_,xm,ym,zm,xn,yn,zn,theta_mp,phi_mp,theta_nm,phi_nm);
g_pm   	=   @(l,l_) exp(-j*k*R_pm(l,l_))./(4*pi*R_pm(l,l_));
Z1_pm   =   @(l,l_) j*k*eta*dot_pm*g_pm(l,l_).*(1-l/lmp).*(1-l_/lnm);
Z2_pm   =   @(l,l_) -j*(eta/k)*(1/(lmp*lnm))*g_pm(l,l_);
Termpm 	=   Quad2(Z1_pm,0,lmp,0,lnm)-Quad2(Z2_pm,0,lmp,0,lnm);
end
%% Term -+
function[Termmp]=Term_mp()
dot_mp  =   sin(theta_mm)*cos(phi_mm)*sin(theta_np)*cos(phi_np)+...
            sin(theta_mm)*sin(phi_mm)*sin(theta_np)*sin(phi_np)+...
            cos(theta_mm)*cos(theta_np);
R_mp    =   @(l,l_) Rmn_mp(l,l_,xm,ym,zm,xn,yn,zn,theta_mm,phi_mm,theta_np,phi_np);
g_mp  	=   @(l,l_) exp(-j*k*R_mp(l,l_))./(4*pi*R_mp(l,l_));
Z1_mp   =   @(l,l_) j*k*eta*dot_mp*g_mp(l,l_).*(1-l/lmm).*(1-l_/lnp);
Z2_mp   =   @(l,l_) -j*(eta/k)*(1/(lmm*lnp))*g_mp(l,l_);
Termmp 	=   Quad2(Z1_mp,0,lmm,0,lnp)-Quad2(Z2_mp,0,lmm,0,lnp);
end
%% R ++
function[R]=Rmn_pp(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn+l*sin(theta_m)*cos(phi_m)-l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn+l*sin(theta_m)*sin(phi_m)-l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn+l*cos(theta_m)-l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% R --
function[R]=Rmn_mm(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn-l*sin(theta_m)*cos(phi_m)+l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn-l*sin(theta_m)*sin(phi_m)+l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn-l*cos(theta_m)+l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% R +-
function[R]=Rmn_pm(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn+l*sin(theta_m)*cos(phi_m)+l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn+l*sin(theta_m)*sin(phi_m)+l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn+l*cos(theta_m)+l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% R -+
function[R]=Rmn_mp(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn-l*sin(theta_m)*cos(phi_m)-l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn-l*sin(theta_m)*sin(phi_m)-l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn-l*cos(theta_m)-l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% Singular term ++ or --
function[I]=SingularTerm(L)
func1  	=   @(l) (1-l/L).*(sqrt(l.^2+a^2)/L+(1-l/L).*0.5.*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l)));
func2  	=   @(l) (1-l/L).*(-sqrt((L-l).^2+a^2)/L+(1-l/L).*0.5.*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l))));
I1    	=   -j*k*L^2/(16*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
func3  	=   @(l) 0.5*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l));
func4  	=   @(l) 0.5*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l)));
I2    	=   -j*k*L^2/(4*pi)+Quad(func3,0,L)/(4*pi)+Quad(func4,0,L)/(4*pi);
I       =   j*k*eta*I1-j*(eta/k)*(1/(L^2))*I2;
end
%% Singular term +- or -+
function[I]=SingularTerm2(L)
func1  	=   @(l) (1/L)*(1-l/L).*(sqrt((L-l).^2+a^2)-sqrt(l.^2+a^2));
func2  	=   @(l) (0.5/L)*(1-l/L).*l.*(log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l))+log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l))));
I1    	=   -j*k*L^2/(16*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
func3  	=   @(l) 0.5*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l));
func4  	=   @(l) 0.5*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l)));
I2    	=   -j*k*L^2/(4*pi)+Quad(func3,0,L)/(4*pi)+Quad(func4,0,L)/(4*pi);
I       =   j*k*eta*I1+j*(eta/k)*(1/(L^2))*I2;
end
end
%%





