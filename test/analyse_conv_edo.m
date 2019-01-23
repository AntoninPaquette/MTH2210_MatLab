%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script effectuant l'analyse de convergence des m�thodes de Runge-Kutta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all


%% D�finition des param�tres pour l'analyse de convergence
nb_eval			=	10;
nb_pas_init		=	100;
nb_pas			=	2.^(0:nb_eval-1) * nb_pas_init;


%% D�finition du probl�me avec sa solution exacte
tspan		=	[0,5];
h			=	(tspan(2)-tspan(1))./(nb_pas);
x0			=	[1,0];
f			=	@my_edo;
sol_exacte	=	@(t) [cos(t*sqrt(10)),-sqrt(10)*sin(t*sqrt(10))];


%% Initialisation des donn�es
erreur_euler_exp	=	nan(nb_eval);
erreur_milieu		=	nan(nb_eval);
erreur_euler_mod	=	nan(nb_eval);
erreur_rk4			=	nan(nb_eval);


%% Appel des diverses m�thodes et calcul des erreurs absolues
for t=1:nb_eval
	[temps , y_euler_exp]		=	crank_nic(f , tspan , x0 , nb_pas(t));
	erreur_euler_exp(t)			=	norm(y_euler_exp - sol_exacte(temps),inf);
% 	[temps , y_milieu]		=	pt_milieu(f , tspan , x0 , nb_pas(t));
% 	erreur_milieu(t)			=	norm(y_milieu - sol_exacte(temps),inf);
% 	[temps , y_euler_mod]	=	euler_mod(f , tspan , x0 , nb_pas(t));
% 	erreur_euler_mod(t)		=	norm(y_euler_mod - sol_exacte(temps),inf);
% 	[temps , y_rk4]				=	rk4(f , tspan , x0 , nb_pas(t));
% 	erreur_rk4(t)					=	norm(y_rk4 - sol_exacte(temps),inf);
end


%% Calcul de l'ordre
ordre_euler_exp	=	log(erreur_euler_exp(1:end-1,:)./erreur_euler_exp(2:end,:))/log(2);
% ordre_milieu	=	log(erreur_milieu(1:end-1,:)./erreur_milieu(2:end,:))/log(2);
% ordre_euler_mod	=	log(erreur_euler_mod(1:end-1,:)./erreur_euler_mod(2:end,:))/log(2);
% ordre_rk4		=	log(erreur_rk4(1:end-1,:)./erreur_rk4(2:end,:))/log(2);


%% Affichage des graphiques de convergence
figure
loglog(h,erreur_euler_exp)
% hold on
% loglog(h,erreur_milieu)
% loglog(h,erreur_euler_mod)
% loglog(h,erreur_rk4)
% xlabel('h')
% ylabel('Erreur absolue')
