close all; clear; clc;
%%
a       =   0.005;
L       =   0.47;   
N       =   31;
%%
tic;
[Data,N,M]=Dipole(N,L);
% PlotDipole(N,M,Data);
I=MM(N,M,Data,a);
toc;
% PlotCurrent(N,Data,I);
%%
if imag(1/I(M))>0
    fprintf('Zin\t=\t%0.2f\t+j\t%0.2f\n',real(1/I(M)),imag(1/I(M)));
else
    fprintf('Zin\t=\t%0.2f\t-j\t%0.2f\n',real(1/I(M)),-imag(1/I(M)));
end
%%
