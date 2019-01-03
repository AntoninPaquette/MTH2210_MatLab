clc
close all
clear

xi=[1;2;4;5];
yi=[1,9,2,11];

x	=	linspace(xi(1),xi(end),200);

Sx =	splinec( xi , yi , x , [3,4] , [-1000,700]);

figure
plot(xi,yi,'or')
hold on
plot(x,Sx,'k')