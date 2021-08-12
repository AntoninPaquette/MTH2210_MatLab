classdef TestNonLineaire < matlab.unittest.TestCase
	
	properties (TestParameter)
		% 2 fonctions avec leurs derivees et racines
		fct		=	{@(x) x^2 - 10, @(x) exp(x) - x^3}
		dfct	=	{@(x) 2*x, @(x) exp(x) - 3*x^2}
		x0		=	{[2,5],[1.5,2.5]}
		
		% 2 fonction avec leur point-fixe, ordre et taux de convergence
		g		=	{@(x) -x^2/10 + x + 1, @(x) -x^2/6 + x + 9/6}
		pt_fixe =	{sqrt(10), 3}
		x0_g	=	{1,1}
		ordre	=	{1,2}
		taux	=	{-2*sqrt(10)/10 + 1, 1/6}
	end
	
	
	methods (Test, ParameterCombination = 'sequential')
		
		function testBissecOrder(testCase,fct,x0)
			% Ordre de convergence pour la methode de la bissection
			
			[~, err] = bissec(fct, x0(1), x0(2), 200 ,1e-12);
			[ordre_app] = order_computation_bissec(err);
			verifyLessThan(testCase,abs(ordre_app-1),0.1);
		end
		
		function testSecanteOrder(testCase,fct,x0)
			% Ordre de convergence pour la methode de la sécante

			[~, err] = secante(fct, x0(1), x0(2), 20 ,1e-12);
			[ordre_app,~] = order_computation_nl(err);
			verifyLessThan(testCase,abs(ordre_app-(1+sqrt(5))/2),0.1);
		end
		
		function testNewton_1DOrder(testCase,fct,dfct,x0)
			% Ordre de convergence pour la methode de Newton

			[~, err] = newton_1D(fct, dfct, x0(1), 20 ,1e-12);
			[ordre_app,~] = order_computation_nl(err,0.4);
			verifyLessThan(testCase,abs(ordre_app-2),0.1);
		end
		
		function testPts_fixes(testCase,g,pt_fixe,x0_g,ordre,taux)
			% Ordre et taux de convergence pour la methode des points-fixes
			
			[app, err] = pts_fixes(g, x0_g, 100 ,1e-12);
			
			verifyLessThan(testCase,abs(app(end)-pt_fixe),1e-11);

			[ordre_app,~] = order_computation_nl(err,0.4);
			verifyLessThan(testCase,abs(ordre_app-ordre),0.1);
			
			[taux_app,~] = taux_computation_nl(err,ordre,0.1);
			verifyLessThan(testCase,abs(taux_app-taux),0.1);
			
		end
	end
	
	methods (Test)
		function testNewton_1DMultiplesRootsOrder(testCase)
			% Ordre et taux de convergence de la methode de Newton pour 
			% racines multiples
			
			fct_mult	=	@(x) x*sin(x)^2;
			dfct_mult	=	@(x) sin(x)^2 + 2*x*sin(x)*cos(x);
			
			[approx_r1 , err_r1] = newton_1D(fct_mult , dfct_mult, 1 , 200 , 1e-12);
			[approx_r2 , err_r2] = newton_1D(fct_mult , dfct_mult, 3 , 200 , 1e-12);
			
			verifyLessThan(testCase,abs(approx_r1(end)),1e-10);
			verifyLessThan(testCase,abs(approx_r2(end)-pi),1e-10);
			
			[ordre_r1,~] = order_computation_nl(err_r1,0.2);
			[ordre_r2,~] = order_computation_nl(err_r2,0.2);
			verifyLessThan(testCase,abs(ordre_r1-1),0.1);
			verifyLessThan(testCase,abs(ordre_r2-1),0.1);
			
			[taux_r1,~] = taux_computation_nl(err_r1,1); % Racine de multiplicite 3
			[taux_r2,~] = taux_computation_nl(err_r2,1); % Racine de multiplicite 2
			verifyLessThan(testCase,abs(taux_r1-(3-1)/3),0.1);
			verifyLessThan(testCase,abs(taux_r2-(2-1)/2),0.1);
		end
		function testNewton_ND_Sans_DerOrder(testCase)
			% Ordre de convergence pour la methode de Newton (avec derivee
			% approximee) pour un systeme non-lineaire
			fct_syst	=	@(x) [5*sin(0.1*x(1)*x(2)) - x(3) , ...
						  x(1)^2 + x(2)^2 + x(3)^2 - 9 , ...
						  x(1) - x(2) - x(3) - 1];
						  
			x0_syst	=	[1;1;1];
			
			[~ , err] = newton_ND_sans_der(fct_syst , x0_syst , 20 , 1e-12);
			
			[ordre_app,~] = order_computation_nl(err,0.4);

			verifyLessThan(testCase,abs(ordre_app-2),0.1);
		end
		function testNewton_ND_Avec_DerOrder(testCase)
			% Ordre de convergence pour la methode de Newton (avec derivee
			% exacte) pour un systeme non-lineaire
			
			fct_syst	=	@(x) [5*sin(0.1*x(1)*x(2)) - x(3) , ...
						  x(1)^2 + x(2)^2 + x(3)^2 - 9 , ...
						  x(1) - x(2) - x(3) - 1];
			
			jac_fct	=	@(x) [5*0.1*x(2)*cos(0.1*x(1)*x(2)) , 5*0.1*x(1)*cos(0.1*x(1)*x(2)) , -1 ;...
							  2*x(1) , 2*x(2) , 2*x(3) ; ...
							  1 , -1 , -1];
						  
			x0_syst	=	[1;1;1];
			
			[~ , err] = newton_ND_avec_der(fct_syst, jac_fct, x0_syst , 20 , 1e-12);
			
			[ordre_app,~] = order_computation_nl(err,0.4);

			verifyLessThan(testCase,abs(ordre_app-2),0.1);
		end
	end
	
end

function [ordre,ordre_app] = order_computation_nl(erreur,varargin)
% Approximation de l'ordre pour les méthodes de resolution de problemes
% non-lineaires

	if nargin>2
		error("Il ne peut y avoir qu'un deuxième argument, la tolérance spécifiée.")
	elseif nargin == 2
		tol = varargin{1};
	else
		tol = 0.2;
	end
	
	ordre_app			=	log(erreur(2:end-1)./erreur(3:end))./log(erreur(1:end-2)./erreur(2:end-1));	
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

function [ordre] = order_computation_bissec(erreur)
% Approximation de l'ordre pour la methode de la bissection avec un
% least-square fit

	% Least-square fit
	nb_iter	=	length(erreur);
	A		=	[ones(nb_iter-2,1) reshape(log(erreur(1:end-2)),[],1)];
	b		=	reshape(log(erreur(2:end-1)),[],1);
	coeff	=	A\b;
	
	ordre	=	coeff(2);

end

function [taux,taux_app] = taux_computation_nl(erreur,ordre,varargin)
% Approximation du taux de convergence

	if nargin>3
		error("Il ne peut y avoir qu'un troisième argument, la tolérance spécifiée.")
	elseif nargin == 3
		tol = varargin{1};
	else
		tol = 0.2;
	end
	
	taux_app			=	erreur(2:end)./(erreur(1:end-1).^ordre);	
	stable_region		=	(taux_app>0) & (abs(gradient(taux_app))<tol) & (abs(del2(taux_app))<2*tol);
	ind_stable_region	=	find(stable_region);
	
	% Sanity check
	if isempty(ind_stable_region)
		error("Il ne pas y avoir de zone asymptotique")
	elseif length(ind_stable_region) < 2
		warning("La zone asymptotique n'est pas très grande")
	elseif any(gradient(ind_stable_region)~=1)
		warning("La zone asymptotique est brisée")
	end
	
	taux = mean(taux_app(ind_stable_region));
end

