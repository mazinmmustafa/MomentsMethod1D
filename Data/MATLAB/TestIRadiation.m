close all; clear; clc;
%%
Ns      =   1E4;
%%
j       =   sqrt(-1);
x       =   logspace(-10,10,Ns);
y       =   (exp(-j*x).*(1+j*x)-1)./x.^2;
%%
figure()
hold on
plot(x,real(y),'-k','LineWidth',2)
y(x<10E-4)  =   0.5;
plot(x,real(y),'-r','LineWidth',1)
hold off
set(gca,'XScale','log')
%%
a       =   1E-1;
b       =   1.4;
ANS     =   f(a,b,'+');
fprintf('(%21.14E %21.14E)\n', real(ANS), imag(ANS));
ANS     =   f(a,b,'-');
fprintf('(%21.14E %21.14E)\n', real(ANS), imag(ANS));
%%
function[F]=f(a,b,s)
j       =   sqrt(-1);
if (abs(a)<1E-4)
    F       =   0.5*exp(j*b);
else
    if s=='+'
        F      =   exp(j*b)*(exp(-j*a).*(1+j*a)-1)./a.^2;
    end
    if s=='-'
        F      =   exp(j*b)*(exp(+j*a).*(1-j*a)-1)./a.^2;
    end
end
end
%%