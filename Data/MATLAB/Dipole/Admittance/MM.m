function[I]=MM(N,M,Data,a)
%%
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
Z       =   zeros(N,N);
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
end
V       =   zeros(N,1);
V(M,1) 	=   1;
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
