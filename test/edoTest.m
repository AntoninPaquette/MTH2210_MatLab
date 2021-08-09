function tests = edoTest
	tests = functiontests(localfunctions);
end

function setup(testCase)
	testCase.TestData.nb_eval		= 12;
	testCase.TestData.nb_pas_init	=	100;
	testCase.TestData.nb_pas	=	2.^(0:testCase.TestData.nb_eval-1) * testCase.TestData.nb_pas_init;
	
	testCase.TestData.tspan	=	[0,5];
	testCase.TestData.h		=	(testCase.TestData.tspan(2)-testCase.TestData.tspan(1))./(testCase.TestData.nb_pas);
	testCase.TestData.x0		=	1;
	testCase.TestData.f1		=	@(t,y) 3;
	testCase.TestData.y1_ex		=	@(t) 3 * t + 1;
	testCase.TestData.f2		=	@(t,y) -4*t;
	testCase.TestData.y2_ex		=	@(t) -2 * t.^2 + 1;
	testCase.TestData.f4		=	@(t,y) 5*t.^3;
	testCase.TestData.y4_ex		=	@(t) 5/4 * t.^4 + 1;
	
	testCase.TestData.f_scalar		=	@(t,y) 2*y -t + 4;
	testCase.TestData.y_scalar_ex	=	@(t) -7/4 +1/2*t + 11/4*exp(2*t);
end

% function teardown(testCase)
% 	
% end

function euler_expExactTest(testCase)
	fct			=	testCase.TestData.f1;
	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	euler_exp(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function euler_modExactTest(testCase)
	fct			=	testCase.TestData.f1;
	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	euler_mod(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	fct			=	testCase.TestData.f2;
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	euler_mod(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);	
end

function pt_milieuExactTest(testCase)
	fct			=	testCase.TestData.f1;
	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	pt_milieu(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	fct			=	testCase.TestData.f2;
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	pt_milieu(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function rk4ExactTest(testCase)
	fct			=	testCase.TestData.f1;
	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	rk4(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	fct			=	testCase.TestData.f2;
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	rk4(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	fct			=	testCase.TestData.f4;
	sol_ex		=	testCase.TestData.y4_ex;
	[temps , y]	=	rk4(fct ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));				
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function euler_expOrdreTest(testCase)
	nb_eval		=	testCase.TestData.nb_eval;
	fct			=	testCase.TestData.f_scalar;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_exp(fct ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-1),tol);
end

function euler_modOrdreTest(testCase)
	nb_eval		=	testCase.TestData.nb_eval;
	fct			=	testCase.TestData.f_scalar;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_mod(fct ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function pt_milieuOrdreTest(testCase)
	nb_eval		=	testCase.TestData.nb_eval;
	fct			=	testCase.TestData.f_scalar;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	pt_milieu(fct ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function rk4OrdreTest(testCase)
	nb_eval		=	testCase.TestData.nb_eval;
	fct			=	testCase.TestData.f_scalar;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	rk4(fct ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-4),tol);
end

function [ordre,ordre_app] = order_computation(erreur,ratio_h,varargin)

	if nargin>3
		error("Il ne peut y avoir qu'un troisième argument, la tolérance spécifiée.")
	elseif nargin == 3
		tol = varargin{1};
	else
		tol = 0.2;
	end
	
	ordre_app			=	log(erreur(1:end-1)./erreur(2:end))/log(ratio_h);	
	stable_region		=	(ordre_app>0) & (abs(gradient(ordre_app))<tol) & (abs(del2(ordre_app))<2*tol);
	ind_stable_region	=	find(stable_region);
	
	% Sanity check
	if isempty(ind_stable_region)
		error("Il ne pas y avoir de zone asymptotique")
	elseif length(ind_stable_region) < 2
		warning("La zone asymptotique n'est pas très grande")
	elseif any(gradient(ind_stable_region)~=1)
		warning("La zone asymptotique est brisée")
	end
	
	ordre = mean(ordre_app(ind_stable_region));
end

%
% clear
% clc
% close all
% 
% tic()
% %% Définition des paramètres pour l'analyse de convergence
% nb_eval			=	12;
% nb_pas_init		=	100;
% nb_pas			=	2.^(0:nb_eval-1) * nb_pas_init;
% 
% 
% %% Définition du problème avec sa solution exacte
% tspan		=	[0,5];
% h			=	(tspan(2)-tspan(1))./(nb_pas);
% x0			=	[1,0];
% f			=	@my_edo;
% sol_exacte	=	@(t) [cos(t*sqrt(10)),-sqrt(10)*sin(t*sqrt(10))];
% 
% 
% %% Initialisation des données
% erreur_euler_exp	=	nan(nb_eval,1);
% erreur_milieu		=	nan(nb_eval,1);
% erreur_euler_mod	=	nan(nb_eval,1);
% erreur_rk4			=	nan(nb_eval,1);
% % erreur_euler_imp	=	nan(nb_eval,1);
% % erreur_crank_nic	=	nan(nb_eval,1);
% 
% 
% 
% %% Appel des diverses méthodes et calcul des erreurs absolues
% for t=1:nb_eval
% 	[temps , y_euler_exp]	=	euler_exp(f , tspan , x0 , nb_pas(t));
% 	erreur_euler_exp(t,:)	=	norm(y_euler_exp - sol_exacte(temps),inf);
% 	[temps , y_milieu]		=	pt_milieu(f , tspan , x0 , nb_pas(t));
% 	erreur_milieu(t,:)		=	norm(y_milieu - sol_exacte(temps),inf);
% 	[temps , y_euler_mod]	=	euler_mod(f , tspan , x0 , nb_pas(t));
% 	erreur_euler_mod(t,:)	=	norm(y_euler_mod - sol_exacte(temps),inf);
% 	[temps , y_rk4]			=	rk4(f , tspan , x0 , nb_pas(t));
% 	erreur_rk4(t,:)			=	norm(y_rk4 - sol_exacte(temps),inf);
% % 	[temps , y_euler_imp]	=	euler_imp(f , tspan , x0 , nb_pas(t));
% % 	erreur_euler_imp(t,:)	=	norm(y_euler_imp - sol_exacte(temps),inf);
% % 	[temps , y_crank_nic]	=	crank_nic(f , tspan , x0 , nb_pas(t));
% % 	erreur_crank_nic(t,:)	=	norm(y_crank_nic - sol_exacte(temps),inf);
% end
% 
% 
% %% Calcul de l'ordre
% ordre_euler_exp	=	log(erreur_euler_exp(1:end-1)./erreur_euler_exp(2:end))/log(2);
% ordre_milieu	=	log(erreur_milieu(1:end-1)./erreur_milieu(2:end))/log(2);
% ordre_euler_mod	=	log(erreur_euler_mod(1:end-1)./erreur_euler_mod(2:end))/log(2);
% ordre_rk4		=	log(erreur_rk4(1:end-1)./erreur_rk4(2:end))/log(2);
% % ordre_euler_imp	=	log(erreur_euler_imp(1:end-1)./erreur_euler_imp(2:end))/log(2);
% % ordre_crank_nic	=	log(erreur_crank_nic(1:end-1)./erreur_crank_nic(2:end))/log(2);
% 
% %% Affichage des graphiques de convergence
% figure
% loglog(h,erreur_euler_exp)
% hold on
% loglog(h,erreur_milieu)
% loglog(h,erreur_euler_mod)
% loglog(h,erreur_rk4)
% % loglog(h,erreur_euler_imp)
% % loglog(h,erreur_crank_nic)
% xlabel('h')
% ylabel('Erreur absolue')
% legend('Euler exp','Pt milieu','Euler mod','Rk4','Euler imp','Crank nic')
% toc()