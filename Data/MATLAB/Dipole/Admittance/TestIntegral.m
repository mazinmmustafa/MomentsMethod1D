close all; clear; clc;
%%
j       =   sqrt(-1);
k       =   2*pi;
L       =   0.05;
a       =   1e-4;
%% Term 1 (new)
fprintf('Term 1\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(alpha-alpha_).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*alpha.*alpha_.*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*alpha.*alpha.*log(((1-alpha)*L+sqrt(a^2+((1-alpha)*L).^2))...
                ./(-alpha*L+sqrt(a^2+(alpha*L).^2)))...
                +alpha.*(sqrt(((1-alpha)*L).^2+a^2)-sqrt((alpha*L).^2+a^2));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%% Term 2 (new)
fprintf('\nTerm 2\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(alpha-alpha_).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*log(((1-alpha)*L+sqrt(a^2+((1-alpha)*L).^2))...
                ./(-alpha*L+sqrt(a^2+(alpha*L).^2)));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi*L^2);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%% Term 3 (new)
fprintf('\nTerm 3\n----------------------------------\n');
R       =   @(alpha,alpha_) sqrt((L^2)*(1-(alpha+alpha_)).^2+a^2);
func1   =   @(alpha,alpha_) -j*k*(L^2)*alpha.*alpha_.*exp(-j*k*R(alpha,alpha_)/2).*sinc(k*R(alpha,alpha_)/2);
func2   =   @(alpha) L*alpha.*(1-alpha).*log((alpha*L+sqrt(a^2+(alpha*L).^2))...
                ./((alpha-1)*L+sqrt(a^2+((alpha-1)*L).^2)))...
                -alpha.*(sqrt(((1-alpha)*L).^2+a^2)-sqrt((alpha*L).^2+a^2));
I1      =   Quad2(func1,0,1,0,1);
I2      =   Quad(func2,0,1);
I       =   (I1+I2)/(4*pi);
fprintf('%0.4e\t+j\t%0.4e (new)\n',real(I),imag(I));
%%
function[y]=sinc(x)
y       =	sin(x)./x;
y(~x)	=	1;
end
%%