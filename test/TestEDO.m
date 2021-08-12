classdef TestEDO < matlab.unittest.TestCase
	
	properties (TestParameter)
		algo		=	{@euler_exp, @euler_imp, @euler_mod, @pt_milieu, @crank_nic, @rk4} 
		algo_degre2	=	{@euler_mod, @pt_milieu, @crank_nic, @rk4}
		algo_degre4 =	{@rk4}
		
		ordre		=	{1,1,2,2,2,4}
	end
	
	methods (Test)
		function testDegre1Exact(testCase,algo)
			% Integration exacte pour solution de degre 1
			
			fct		=	@(t,y) 3;
			y0		=	1;
			y_ex	=	@(t) 3 * t + 1;
			tspan	=	[0,5];
			[temps , y]	=	algo(fct ,tspan, y0, 100);
			erreur_rel	=	norm(y - y_ex(temps),inf)/norm(y_ex(temps),inf);
			verifyLessThan(testCase,erreur_rel,1e-14);
		end
		
		function testDegre2Exact(testCase,algo_degre2)
			% Integration exacte pour solution de degre 2

			fct		=	@(t,y) -4*t;
			y0		=	1;
			y_ex	=	@(t) -2 * t.^2 + 1;
			tspan	=	[0,5];
			[temps , y]	=	algo_degre2(fct ,tspan, y0, 100);
			erreur_rel	=	norm(y - y_ex(temps),inf)/norm(y_ex(temps),inf);
			
			verifyLessThan(testCase,erreur_rel,1e-14);
		end
		
		function testDegre4Exact(testCase,algo_degre4)
			% Integration exacte pour solution de degre 4

			fct		=	@(t,y) 5*t.^3;
			y0		=	1;
			y_ex	=	@(t) 5/4 * t.^4 + 1;
			tspan	=	[0,5];
			[temps , y]	=	algo_degre4(fct ,tspan, y0, 100);
			erreur_rel	=	norm(y - y_ex(temps),inf)/norm(y_ex(temps),inf);
			
			verifyLessThan(testCase,erreur_rel,1e-14);
		end
	end
	
	methods (Test, ParameterCombination = 'sequential')
		function testOrdreScalar(testCase,algo,ordre)
			% Ordre de convergence pour une EDO (pour toutes les méthodes)
			
			fct		=	@(t,y) 2*y -t + 4;
			y0		=	1;
			y_ex	=	@(t) -7/4 +1/2*t + 11/4*exp(2*t);
			tspan	=	[0,5];
			
			nb_eval		=	8;
			nb_pas_init	=	100;
			nb_pas		=	2.^(0:nb_eval-1) * nb_pas_init;
			erreur		=	nan(nb_eval,1);
			
			for t=1:nb_eval
				[temps , y]	=	algo(fct ,tspan, y0, nb_pas(t));
				erreur(t)	=	norm(y - y_ex(temps),inf);
			end
			
			tol	=	0.2;
			[ordre_app,~] = order_computation(erreur,2,tol);
			
			verifyLessThan(testCase,abs(ordre_app-ordre),tol);
		end
		
		function testOrdreSystem(testCase,algo,ordre)
			% Ordre de convergence pour un systeme EDO (pour toutes les methodes)

			fct		=	@(t,y) [-2,1;1,-2]*y(:) + [2*exp(-t);3*t];
			y0		=	[2;3];
			y_ex	=	@(t) -7/6*[1;-1]*exp(-3*t) + 4*[1;1]*exp(-t) + ...
							1/2*[1;-1]*exp(-t) + [1;1]*t.*exp(-t) + [1;2]*t -1/3*[4;5];
			tspan	=	[0,5];
			
			nb_eval		=	8;
			nb_pas_init	=	100;
			nb_pas		=	2.^(0:nb_eval-1) * nb_pas_init;
			erreur		=	nan(nb_eval,1);
			
			for t=1:nb_eval
				[temps , y]	=	algo(fct ,tspan, y0, nb_pas(t));
				erreur(t)	=	norm(y - y_ex(temps),1);
			end
			
			tol	=	0.2;
			[ordre_app,~] = order_computation(erreur,2,tol)
			
			verifyLessThan(testCase,abs(ordre_app-ordre),tol);
		end
	end
end


function [ordre,ordre_app] = order_computation(erreur,ratio_h,varargin)
% Approximation de l'ordre de convergence

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


