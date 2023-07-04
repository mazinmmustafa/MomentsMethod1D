close all; clear; clc;
%%

Data = load('data1.dat');


figure()
plot(Data(:,1), Data(:,2:end))
##plot(Data(:,1), Data(:,2:end),'o')


%%