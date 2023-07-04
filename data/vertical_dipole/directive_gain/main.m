close all; clear; clc;

Data = load("data1.dat");

polar(Data(:,1),Data(:,3))