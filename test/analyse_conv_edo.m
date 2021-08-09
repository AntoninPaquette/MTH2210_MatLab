%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script effectuant l'analyse de convergence des méthodes de Runge-Kutta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

tic()
%% Définition des paramètres pour l'analyse de convergence
nb_eval			=	12;
nb_pas_init		=	100;
nb_pas			=	2.^(0:nb_eval-1) * nb_pas_init;


%% Définition du problème avec sa solution exacte
tspan		=	[0,5];
h			=	(tspan(2)-tspan(1))./(nb_pas);
x0			=	[1,0];
f			=	@my_edo;
sol_exacte	=	@(t) [cos(t*sqrt(10)),-sqrt(10)*sin(t*sqrt(10))];


%% Initialisation des données
erreur_euler_exp	=	nan(nb_eval,1);
erreur_milieu		=	nan(nb_eval,1);
erreur_euler_mod	=	nan(nb_eval,1);
erreur_rk4			=	nan(nb_eval,1);
% erreur_euler_imp	=	nan(nb_eval,1);
% erreur_crank_nic	=	nan(nb_eval,1);



%% Appel des diverses méthodes et calcul des erreurs absolues
for t=1:nb_eval
	[temps , y_euler_exp]	=	euler_exp(f , tspan , x0 , nb_pas(t));
	erreur_euler_exp(t,:)	=	norm(y_euler_exp - sol_exacte(temps),inf);
	[temps , y_milieu]		=	pt_milieu(f , tspan , x0 , nb_pas(t));
	erreur_milieu(t,:)		=	norm(y_milieu - sol_exacte(temps),inf);
	[temps , y_euler_mod]	=	euler_mod(f , tspan , x0 , nb_pas(t));
	erreur_euler_mod(t,:)	=	norm(y_euler_mod - sol_exacte(temps),inf);
	[temps , y_rk4]			=	rk4(f , tspan , x0 , nb_pas(t));
	erreur_rk4(t,:)			=	norm(y_rk4 - sol_exacte(temps),inf);
% 	[temps , y_euler_imp]	=	euler_imp(f , tspan , x0 , nb_pas(t));
% 	erreur_euler_imp(t,:)	=	norm(y_euler_imp - sol_exacte(temps),inf);
% 	[temps , y_crank_nic]	=	crank_nic(f , tspan , x0 , nb_pas(t));
% 	erreur_crank_nic(t,:)	=	norm(y_crank_nic - sol_exacte(temps),inf);
end


%% Calcul de l'ordre
ordre_euler_exp	=	log(erreur_euler_exp(1:end-1)./erreur_euler_exp(2:end))/log(2);
ordre_milieu	=	log(erreur_milieu(1:end-1)./erreur_milieu(2:end))/log(2);
ordre_euler_mod	=	log(erreur_euler_mod(1:end-1)./erreur_euler_mod(2:end))/log(2);
ordre_rk4		=	log(erreur_rk4(1:end-1)./erreur_rk4(2:end))/log(2);
% ordre_euler_imp	=	log(erreur_euler_imp(1:end-1)./erreur_euler_imp(2:end))/log(2);
% ordre_crank_nic	=	log(erreur_crank_nic(1:end-1)./erreur_crank_nic(2:end))/log(2);

%% Affichage des graphiques de convergence
figure
loglog(h,erreur_euler_exp)
hold on
loglog(h,erreur_milieu)
loglog(h,erreur_euler_mod)
loglog(h,erreur_rk4)
% loglog(h,erreur_euler_imp)
% loglog(h,erreur_crank_nic)
xlabel('h')
ylabel('Erreur absolue')
legend('Euler exp','Pt milieu','Euler mod','Rk4','Euler imp','Crank nic')
toc()