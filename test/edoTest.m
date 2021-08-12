function tests = edoTest
	tests = functiontests(localfunctions);
end

function setup(testCase)
	testCase.TestData.nb_eval		= 8;
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
	
	testCase.TestData.x0_syst		=	[2;3];
	testCase.TestData.f_syst		=	@(t,y) [-2,1;1,-2]*y(:) + [2*exp(-t);3*t];
	testCase.TestData.y_syst_ex		=	@(t) -7/6*[1;-1]*exp(-3*t) + 4*[1;1]*exp(-t) + ...
											 1/2*[1;-1]*exp(-t) + [1;1]*t.*exp(-t) + ...
											 [1;2]*t -1/3*[4;5];
end

% function teardown(testCase)
% 	
% end

function euler_expExactTest(testCase)

	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	euler_exp(testCase.TestData.f1 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function euler_impExactTest(testCase)

	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	euler_imp(testCase.TestData.f1 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function euler_modExactTest(testCase)

	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	euler_mod(testCase.TestData.f1 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	euler_mod(testCase.TestData.f2 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);	
end

function pt_milieuExactTest(testCase)

	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	pt_milieu(testCase.TestData.f1 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	pt_milieu(testCase.TestData.f2 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function crank_nicExactTest(testCase)

	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	crank_nic(testCase.TestData.f1 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	pt_milieu(testCase.TestData.f2 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function rk4ExactTest(testCase)

	sol_ex		=	testCase.TestData.y1_ex;
	[temps , y]	=	rk4(testCase.TestData.f1 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	sol_ex		=	testCase.TestData.y2_ex;
	[temps , y]	=	rk4(testCase.TestData.f2 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));					
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
	
	sol_ex		=	testCase.TestData.y4_ex;
	[temps , y]	=	rk4(testCase.TestData.f4 ,testCase.TestData.tspan , ...
							testCase.TestData.x0 , testCase.TestData.nb_pas(1));				
	erreur_rel	=	norm(y - sol_ex(temps),inf)/norm(sol_ex(temps),inf);
	verifyLessThan(testCase,erreur_rel,1e-14);
end

function euler_expOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_exp(testCase.TestData.f_scalar ,testCase.TestData.tspan , ...
								  testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-1),tol);
end

function euler_impOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_imp(testCase.TestData.f_scalar ,testCase.TestData.tspan , ...
								  testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-1),tol);
end

function euler_modOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_mod(testCase.TestData.f_scalar ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function pt_milieuOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	pt_milieu(testCase.TestData.f_scalar ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function crank_nicOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	crank_nic(testCase.TestData.f_scalar ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function rk4OrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_scalar_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	rk4(testCase.TestData.f_scalar ,testCase.TestData.tspan , ...
								testCase.TestData.x0 , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),inf);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-4),tol);
end

function euler_expSystOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_syst_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_exp(testCase.TestData.f_syst ,testCase.TestData.tspan , ...
								  testCase.TestData.x0_syst , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),1);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-1),tol);
end

function euler_impSystOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_syst_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_imp(testCase.TestData.f_syst ,testCase.TestData.tspan , ...
								  testCase.TestData.x0_syst , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),1);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-1),tol);
end

function euler_modSystOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_syst_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	euler_mod(testCase.TestData.f_syst ,testCase.TestData.tspan , ...
								testCase.TestData.x0_syst , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),1);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function pt_milieuSystOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_syst_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	pt_milieu(testCase.TestData.f_syst ,testCase.TestData.tspan , ...
								testCase.TestData.x0_syst , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),1);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function crank_nicSystOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_syst_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	crank_nic(testCase.TestData.f_syst ,testCase.TestData.tspan , ...
								  testCase.TestData.x0_syst , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),1);
	end
	
	tol	=	0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-2),tol);
end

function rk4SystOrdreTest(testCase)

	nb_eval		=	testCase.TestData.nb_eval;
	sol_ex		=	testCase.TestData.y_syst_ex;
	erreur		=	nan(nb_eval,1);
	
	for t=1:nb_eval
		[temps , y]	=	rk4(testCase.TestData.f_syst ,testCase.TestData.tspan , ...
								testCase.TestData.x0_syst , testCase.TestData.nb_pas(t));
		erreur(t)	=	norm(y - sol_ex(temps),1);
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
	
	ordre_app			=	log(erreur(1:end-1)./erreur(2:end))./log(ratio_h);	
	stable_region		=	(ordre_app>0) & (abs(gradient(ordre_app))<tol) & (abs(del2(ordre_app))<2*tol);
	ind_stable_region	=	find(stable_region);
	
	% Sanity check
	if isempty(ind_stable_region)
		error("Il n'y a pas de zone asymptotique")
	elseif length(ind_stable_region) < 2
		warning("La zone asymptotique n'est pas très grande")
	elseif any(gradient(ind_stable_region)~=1)
		warning("La zone asymptotique est brisée")
	end
	
	ordre = mean(ordre_app(ind_stable_region));
end
