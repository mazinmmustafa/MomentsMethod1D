close all; clear; clc;
%%
j       =   sqrt(-1);
k       =   2*pi;
L       =   0.05;
a       =   1e-4;
%% Term 1 (new)
fprintf('Tern 1\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(alpha-alpha_).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*alpha.*alpha_.*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*alpha.*alpha.*log(((1-alpha)*L+sqrt(a^2+((1-alpha)*L).^2))...
                ./(-alpha*L+sqrt(a^2+(alpha*L).^2)))...
                +alpha.*(sqrt(((1-alpha)*L).^2+a^2)-sqrt((alpha*L).^2+a^2));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%% Term 1 (old)
func1  	=   @(l) (1/L)*(1-l/L).*(sqrt(l.^2+a^2)-sqrt((L-l).^2+a^2));
func2  	=   @(l) 0.5*((1-l/L).^2).*(log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l))+log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l))));
I    	=   -j*k*L^2/(16*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (old)\n',real(I),imag(I))
R       =   @(l,l_) sqrt((l-l_).^2+a^2);
func  	=   @(l,l_) (1-l/L).*(1-l_/L).*exp(-j*k*R(l,l_))./R(l,l_);
I    	=   integral2(func,0,L,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (MATLAB)\n',real(I),imag(I))
fprintf('----------------------------------\n');
%% Term 2 (new)
fprintf('\nTern 2\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(alpha-alpha_).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*log(((1-alpha)*L+sqrt(a^2+((1-alpha)*L).^2))...
                ./(-alpha*L+sqrt(a^2+(alpha*L).^2)));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%% Term 2 (old)
func1  	=   @(l) 0.5*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l));
func2  	=   @(l) 0.5*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l)));
I    	=   -j*k*L^2/(4*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (old)\n',real(I),imag(I))
R       =   @(l,l_) sqrt((l-l_).^2+a^2);
func  	=   @(l,l_) exp(-j*k*R(l,l_))./R(l,l_);
I    	=   integral2(func,0,L,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (MATLAB)\n',real(I),imag(I))
fprintf('----------------------------------\n');
%% Term 3 (new)
fprintf('\nTern 3\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(1-(alpha+alpha_)).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*alpha.*alpha_.*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*alpha.*(1-alpha).*log((alpha*L+sqrt(a^2+(alpha*L).^2))...
                ./((alpha-1)*L+sqrt(a^2+((alpha-1)*L).^2)))...
                -alpha.*(sqrt(((1-alpha)*L).^2+a^2)-sqrt((alpha*L).^2+a^2));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%% Term 3 (old)
func1  	=   @(l) (1/L)*(1-l/L).*(sqrt((L-l).^2+a^2)-sqrt(l.^2+a^2));
func2  	=   @(l) (0.5/L)*(1-l/L).*l.*(log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l))+log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l))));
I    	=   -j*k*L^2/(16*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (old)\n',real(I),imag(I))
R       =   @(l,l_) sqrt((L-(l+l_)).^2+a^2);
func  	=   @(l,l_) (1-l/L).*(1-l_/L).*exp(-j*k*R(l,l_))./R(l,l_);
I    	=   integral2(func,0,L,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (MATLAB)\n',real(I),imag(I))
fprintf('----------------------------------\n');
%% Term 4 (new)
fprintf('\nTern 4\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(1-(alpha+alpha_)).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*log((alpha*L+sqrt(a^2+(alpha*L).^2))...
                ./((alpha-1)*L+sqrt(a^2+((alpha-1)*L).^2)));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%% Term 4 (old)
func1  	=   @(l) 0.5*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l));
func2  	=   @(l) 0.5*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l)));
I    	=   -j*k*L^2/(4*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (old)\n',real(I),imag(I))
R       =   @(l,l_) sqrt((L-(l+l_)).^2+a^2);
func  	=   @(l,l_) exp(-j*k*R(l,l_))./R(l,l_);
I    	=   integral2(func,0,L,0,L)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (MATLAB)\n',real(I),imag(I))
fprintf('----------------------------------\n');
%%
function[y]=sinc(x)
y       =	sin(x)./x;
y(~x)	=	1;
end
%%