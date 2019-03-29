%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script effectuant l'analyse de convergence des méthodes de résolution de
% problèmes non-linéaires
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% Vérification de la fonction lagrange pour polynôme 1

fct1	=	@(x) (x-1).*(x+2.5).*(x-pi).*(x+11).^2;

xi		=	[-10,-2,4,pi,-exp(1),0];
yi		=	fct1(xi);
xfin	=	linspace(min(xi),max(xi),1000);

y_inter		=	lagrange(xi,yi,xfin);
y_exacte	=	fct1(xfin);

err_rel1	=	norm(y_exacte - y_inter)/norm(y_exacte);

%% Vérification de la fonction lagrange pour polynôme 2

fct2	=	@(x) (x-2).^2 .* (x+5).^3 .* (x-exp(1)).^4;

xi		=	linspace(-2,2,10);
yi		=	fct2(xi);
xfin	=	linspace(min(xi),max(xi),1000);

y_inter		=	lagrange(xi,yi,xfin);
y_exacte	=	fct2(xfin);

err_rel2	=	norm(y_exacte - y_inter)/norm(y_exacte);


%% Vérification avec un fonction quelconque

fct3	=	@(x) cos(x);

degre	=	2;
nb_pts	=	degre + 1;
nb_loop		=	10;

x_interet	=	1/2^nb_loop;
erreur		=	nan(nb_loop,1);

for t=1:nb_loop
	xi	=	linspace(-1/2^(t-1),1/2^(t-1),nb_pts);
	yi	=	fct3(xi);
	
	y_inter		=	lagrange(xi,yi,x_interet);
	erreur(t)	=	abs(fct3(x_interet) - y_inter);
end
	
ordre	=	log(erreur(1:end-1)./erreur(2:end))/log(2);


%% Vérification de la fonction splinec

fct4		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
d_fct4		=	@(x) -12*x.^2 + 2*pi*x + 11;
d2_fct4		=	@(x) -24*x + 2*pi;

a		=	-5;
b		=	4;
nb_pts	=	3;

xi	=	linspace(a,b,nb_pts);
yi	=	fct4(xi);
x	=	linspace(a,b,1000);
y_exacte	=	fct4(x);

[ Sx ] = splinec( xi , yi , x , [4,4] , [d_fct4(a),d_fct4(b)]);

err_rel4	=	norm(y_exacte - Sx)/norm(y_exacte);


