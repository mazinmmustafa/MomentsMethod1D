function[I]=MM(N,Data,a,E_TM,E_TE,theta_i,phi_i)
%%
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
Z       =   zeros(N,N);
V       =   zeros(N,1);
for m=1:N
    xm      =   Data(m,1);   
    ym      =   Data(m,2);
    zm      =   Data(m,3);
    xmm     =   Data(m,4);   
    ymm  	=   Data(m,5);
    zmm  	=   Data(m,6);
    xmp     =   Data(m,7);   
    ymp     =   Data(m,8);   
    zmp  	=   Data(m,9);
    Lmm     =   sqrt((xm-xmm)^2+(ym-ymm)^2+(zm-zmm)^2);
    Lmp     =   sqrt((xm-xmp)^2+(ym-ymp)^2+(zm-zmp)^2);
    for n=1:N
        xn      =   Data(n,1);   
        yn      =   Data(n,2);
        zn      =   Data(n,3);
        xnm     =   Data(n,4);   
        ynm     =   Data(n,5);   
        znm  	=   Data(n,6);
        xnp     =   Data(n,7);   
        ynp     =   Data(n,8);   
        znp  	=   Data(n,9);
        if m==n
            Z(m,n)  =   j*k*eta*(I1(Lmp)+Z1pm+Z1mp+I1(Lmm))-j*(eta/k)*(I2(Lmp)-Z2pm-Z2mp+I2(Lmm));
        elseif m==n+1
            Z(m,n)  =   j*k*eta*(Z1pp+Z1pm+I3(Lmm)+Z1mm)-j*(eta/k)*(Z2pp-Z2pm-I2(Lmm)+Z2mm);
        elseif n==m+1
            Z(m,n)  =   j*k*eta*(Z1pp+I3(Lmp)+Z1mp+Z1mm)-j*(eta/k)*(Z2pp-I2(Lmp)-Z2mp+Z2mm);
        else
            Z(m,n)  =   j*k*eta*(Z1pp+Z1pm+Z1mp+Z1mm)-j*(eta/k)*(Z2pp-Z2pm-Z2mp+Z2mm);
        end
    end
    factor_p    =   exp(j*k*(xmp*sin(theta_i)*cos(phi_i)+ymp*sin(theta_i)*sin(phi_i)+zmp*cos(theta_i)));
    factor_m    =   exp(j*k*(xmm*sin(theta_i)*cos(phi_i)+ymm*sin(theta_i)*sin(phi_i)+zmm*cos(theta_i)));
    chi_p   =   @(alpha) alpha.*exp(-j*alpha*k*((xmp-xm)*sin(theta_i)*cos(phi_i)+(ymp-ym)*sin(theta_i)*sin(phi_i)+(zmp-zm)*cos(theta_i)));
    chi_m   =   @(alpha) alpha.*exp(+j*alpha*k*((xm-xmm)*sin(theta_i)*cos(phi_i)+(ym-ymm)*sin(theta_i)*sin(phi_i)+(zm-zmm)*cos(theta_i)));
    term_p  =   factor_p*Quad(@(alpha)chi_p(alpha),0,1);
    term_m  =   factor_m*Quad(@(alpha)chi_m(alpha),0,1);
    dot1_p 	=   (xmp-xm)*cos(theta_i)*cos(phi_i)+(ymp-ym)*cos(theta_i)*sin(phi_i)-(zmp-zm)*sin(theta_i);
    dot1_m 	=   (xm-xmm)*cos(theta_i)*cos(phi_i)+(ym-ymm)*cos(theta_i)*sin(phi_i)-(zm-zmm)*sin(theta_i);
    dot2_p  =   -(xmp-xm)*sin(phi_i)+(ymp-ym)*cos(phi_i);
    dot2_m  =   -(xm-xmm)*sin(phi_i)+(ym-ymm)*cos(phi_i);
    V(m)    =  E_TM*(dot1_p*term_p+dot1_m*term_m)+...
               E_TE*(dot2_p*term_p+dot2_m*term_m);
end
I       =   Z\V;
%% Term ++
function[I]=Z1pp()
Rx     	=   @(alpha,alpha_) xmp-xnp-alpha*(xmp-xm)+alpha_*(xnp-xn);
Ry     	=   @(alpha,alpha_) ymp-ynp-alpha*(ymp-ym)+alpha_*(ynp-yn);
Rz     	=   @(alpha,alpha_) zmp-znp-alpha*(zmp-zm)+alpha_*(znp-zn);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) alpha.*alpha_.*g(alpha,alpha_);
dot     =   (xmp-xm)*(xnp-xn)+(ymp-ym)*(ynp-yn)+(zmp-zm)*(znp-zn);
I     	=   dot*Quad2(func,0,1,0,1);
end
function[I]=Z2pp()
Rx     	=   @(alpha,alpha_) xmp-xnp-alpha*(xmp-xm)+alpha_*(xnp-xn);
Ry     	=   @(alpha,alpha_) ymp-ynp-alpha*(ymp-ym)+alpha_*(ynp-yn);
Rz     	=   @(alpha,alpha_) zmp-znp-alpha*(zmp-zm)+alpha_*(znp-zn);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) g(alpha,alpha_);
I     	=   Quad2(func,0,1,0,1);
end
%% Term +-
function[I]=Z1pm()
Rx     	=   @(alpha,alpha_) xmp-xnm-alpha*(xmp-xm)-alpha_*(xn-xnm);
Ry     	=   @(alpha,alpha_) ymp-ynm-alpha*(ymp-ym)-alpha_*(yn-ynm);
Rz     	=   @(alpha,alpha_) zmp-znm-alpha*(zmp-zm)-alpha_*(zn-znm);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) alpha.*alpha_.*g(alpha,alpha_);
dot     =   (xmp-xm)*(xn-xnm)+(ymp-ym)*(yn-ynm)+(zmp-zm)*(zn-znm);
I     	=   dot*Quad2(func,0,1,0,1);
end
function[I]=Z2pm()
Rx     	=   @(alpha,alpha_) xmp-xnm-alpha*(xmp-xm)-alpha_*(xn-xnm);
Ry     	=   @(alpha,alpha_) ymp-ynm-alpha*(ymp-ym)-alpha_*(yn-ynm);
Rz     	=   @(alpha,alpha_) zmp-znm-alpha*(zmp-zm)-alpha_*(zn-znm);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) g(alpha,alpha_);
I     	=   Quad2(func,0,1,0,1);
end
%% Term -+
function[I]=Z1mp()
Rx     	=   @(alpha,alpha_) xmm-xnp+alpha*(xm-xmm)+alpha_*(xn-xnm);
Ry     	=   @(alpha,alpha_) ymm-ynp+alpha*(ym-ymm)+alpha_*(yn-ynm);
Rz     	=   @(alpha,alpha_) zmm-znp+alpha*(zm-zmm)+alpha_*(zn-znm);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) alpha.*alpha_.*g(alpha,alpha_);
dot     =   (xm-xmm)*(xnp-xn)+(ym-ymm)*(ynp-yn)+(zm-zmm)*(znp-zn);
I    	=   dot*Quad2(func,0,1,0,1);
end
function[I]=Z2mp()
Rx     	=   @(alpha,alpha_) xmm-xnp+alpha*(xm-xmm)+alpha_*(xn-xnm);
Ry     	=   @(alpha,alpha_) ymm-ynp+alpha*(ym-ymm)+alpha_*(yn-ynm);
Rz     	=   @(alpha,alpha_) zmm-znp+alpha*(zm-zmm)+alpha_*(zn-znm);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) g(alpha,alpha_);
I     	=   Quad2(func,0,1,0,1);
end
%% Term --
function[I]=Z1mm()
Rx     	=   @(alpha,alpha_) xmm-xnm+alpha*(xm-xmm)-alpha_*(xn-xnm);
Ry     	=   @(alpha,alpha_) ymm-ynm+alpha*(ym-ymm)-alpha_*(yn-ynm);
Rz     	=   @(alpha,alpha_) zmm-znm+alpha*(zm-zmm)-alpha_*(zn-znm);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) alpha.*alpha_.*g(alpha,alpha_);
dot     =   (xm-xmm)*(xn-xnm)+(ym-ymm)*(yn-ynm)+(zm-zmm)*(zn-znm);
I   	=   dot*Quad2(func,0,1,0,1);
end
function[I]=Z2mm()
Rx     	=   @(alpha,alpha_) xmm-xnm+alpha*(xm-xmm)-alpha_*(xn-xnm);
Ry     	=   @(alpha,alpha_) ymm-ynm+alpha*(ym-ymm)-alpha_*(yn-ynm);
Rz     	=   @(alpha,alpha_) zmm-znm+alpha*(zm-zmm)-alpha_*(zn-znm);
R       =   @(alpha,alpha_) sqrt(Rx(alpha,alpha_).^2+Ry(alpha,alpha_).^2+Rz(alpha,alpha_).^2+a^2);
g       =   @(alpha,alpha_) exp(-j*k*R(alpha,alpha_))./(4*pi*R(alpha,alpha_));
func    =   @(alpha,alpha_) g(alpha,alpha_);
I     	=   Quad2(func,0,1,0,1);
end
%%
function[I]=I1(L)
R       =   @(alpha,alpha_) sqrt((L^2)*(alpha-alpha_).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*alpha.*alpha_.*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*alpha.*alpha.*log(((1-alpha)*L+sqrt(a^2+((1-alpha)*L).^2))...
                ./(-alpha*L+sqrt(a^2+(alpha*L).^2)))...
                +alpha.*(sqrt(((1-alpha)*L).^2+a^2)-sqrt((alpha*L).^2+a^2));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
end
%%
function[I]=I2(L)
R       =   @(alpha,alpha_) sqrt((L^2)*(alpha-alpha_).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*log(((1-alpha)*L+sqrt(a^2+((1-alpha)*L).^2))...
                ./(-alpha*L+sqrt(a^2+(alpha*L).^2)));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi*L^2);
end
%%
function[I]=I3(L)
R       =   @(alpha,alpha_) sqrt((L^2)*(1-(alpha+alpha_)).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*alpha.*alpha_.*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*alpha.*(1-alpha).*log((alpha*L+sqrt(a^2+(alpha*L).^2))...
                ./((alpha-1)*L+sqrt(a^2+((alpha-1)*L).^2)))...
                -alpha.*(sqrt(((1-alpha)*L).^2+a^2)-sqrt((alpha*L).^2+a^2));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
end
%%
end
