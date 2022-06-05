close all; clear; clc;
%% 1D
tic;
Quad(@(x)sin(sqrt(x).*cos(x)),0,pi)
toc;
integral(@(x)sin(sqrt(x).*cos(x)),0,pi)
%% 2D
tic;
Quad2(@(x,y)sin(sqrt(x).*cos(y)),0,2*pi,0,pi/2)
toc;
integral2(@(x,y)sin(sqrt(x).*cos(y)),0,2*pi,0,pi/2)
%%