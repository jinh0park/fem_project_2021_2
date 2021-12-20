clear all;
close all;
clc;
%% 
x=0:0.01:1;
y1=0.7*x.^3/12 + (-3*0.7+0.7^3)*x.^2/48;
y2=-(x-0.15).^4/24;
y3=(x-0.85).^4/24;

y=y1;
y(16:end) = y(16:end)+y2(16:end);
y(86:end) = y(86:end)+y3(86:end);

plot(x,y)