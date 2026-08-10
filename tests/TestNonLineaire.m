classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture( ...
        '../source')}) TestNonLineaire < matlab.unittest.TestCase
	
	properties (TestParameter)
		% 2 fonctions avec leurs derivees et racines
		fct		=	{@(x) x^2 - 10, @(x) exp(x) - x^3 - (exp(pi) - pi^3)}
		dfct	=	{@(x) 2*x, @(x) exp(x) - 3*x^2}
		x0		=	{[2,4],[2.75,3.25]}
		racine	=	{sqrt(10),pi}
		
		% 2 fonction avec leur point-fixe, ordre et taux de convergence
		g		=	{@(x) -x^2/10 + x + 1, @(x) -x^2/6 + x + 9/6}
		pt_fixe =	{sqrt(10), 3}
		x0_g	=	{1,1}
		ordre	=	{1,2}
		taux	=	{-2*sqrt(10)/10 + 1, 1/6}
		
		% Verification pour size des inputs
		f_size	=	{@(x) [5*sin(0.1*x(1)*x(2)) - x(3) - (5*sin(-0.2)-5); ...
						  x(1)^2 + x(2)^2 + x(3)^2 - 30; ...
						  x(1) - x(2) - x(3) + 2],...
					 @(x) [5*sin(0.1*x(1)*x(2)) - x(3) - (5*sin(-0.2)-5); ...
						  x(1)^2 + x(2)^2 + x(3)^2 - 30; ...
						  x(1) - x(2) - x(3) + 2]'}
		x0_size	=	{[1.5,-3,4],[1.5,-3,4]'}
	end
	
	
	methods (Test, ParameterCombination = 'sequential')
		
		function testBissec(testCase,fct,x0,racine)
			% Ordre de convergence pour la methode de la bissection
			
			[app, err] = bissec(fct, x0(1), x0(2), 200 ,1e-12);
			
			verifyLessThan(testCase,abs(app(end)-racine),1e-11);
			
			[ordre_app] = order_computation_bissec(err);
			verifyLessThan(testCase,abs(ordre_app-1),0.1);
		end
		
		function testSecante(testCase,fct,x0,racine)
			% Ordre de convergence pour la methode de la secante

			[app, err] = secante(fct, x0(1), x0(2), 20 ,1e-12);
			
			verifyLessThan(testCase,abs(app(end)-racine),1e-11);
			
			[ordre_app] = order_computation_nl(err,0.4);
			verifyLessThan(testCase,abs(ordre_app-(1+sqrt(5))/2),0.1);
		end
		
		function testNewton_1D(testCase,fct,dfct,x0,racine)
			% Ordre de convergence pour la methode de Newton

			[app, err] = newton_1D(fct, dfct, x0(1), 20 ,1e-12);
			
			verifyLessThan(testCase,abs(app(end)-racine),1e-11);
			
			[ordre_app] = order_computation_nl(err,0.4);
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
			fct_syst	=	@(x) [5*sin(0.1*x(1)*x(2)) - x(3) - (5*sin(-0.2)-5); ...
						  x(1)^2 + x(2)^2 + x(3)^2 - 30; ...
						  x(1) - x(2) - x(3) + 2]';
			
			racine_syst		=	[1;-2;5];		  
			x0_syst	=	[1.5;-3;4];
			
			[app , err] = newton_ND_sans_der(fct_syst, 1e-3, x0_syst , 20 , 1e-12);
			
			verifyLessThan(testCase,norm(app(:,end)-reshape(racine_syst,[],1)),1e-11);
			
			[ordre_app,~] = order_computation_nl(err,0.3);
			verifyLessThan(testCase,abs(ordre_app-2),0.1);
		end
		function testNewton_ND_Avec_DerOrder(testCase)
			% Ordre de convergence pour la methode de Newton (avec derivee
			% exacte) pour un systeme non-lineaire
			
			fct_syst	=	@(x) [5*sin(0.1*x(1)*x(2)) - x(3) - (5*sin(-0.2)-5); ...
						  x(1)^2 + x(2)^2 + x(3)^2 - 30; ...
						  x(1) - x(2) - x(3) + 2];
			
			jac_fct	=	@(x) [5*0.1*x(2)*cos(0.1*x(1)*x(2)) , 5*0.1*x(1)*cos(0.1*x(1)*x(2)) , -1 ;...
							  2*x(1) , 2*x(2) , 2*x(3) ; ...
							  1 , -1 , -1];
						  
			racine_syst		=	[1;-2;5];		  
			x0_syst	=	[1.5,-3,4];
			
			[app , err] = newton_ND_avec_der(fct_syst, jac_fct, x0_syst , 20 , 1e-12);
			
			verifyLessThan(testCase,norm(app(:,end)-reshape(racine_syst,[],1)),1e-11);
			
			[ordre_app,~] = order_computation_nl(err,0.3);
			verifyLessThan(testCase,abs(ordre_app-2),0.1);
		end
		
		function testSizeNewton_ND_Sans_Der(testCase,f_size,x0_size)
			% Robustesse face a l'orientation des vecteurs (rangee ou 
			% colonne) pour les arguments d'entrees (sauf pour la matrice jacobienne) 
					
			[app , err] = newton_ND_sans_der(f_size, 1e-3, x0_size, 20, 1e-12);
			
			nb_iter		=	length(err);
			
			verifySize(testCase,app,[3,nb_iter]);
			verifySize(testCase,err,[1,nb_iter]);
			
		end
		
		function testSizeNewton_ND_Avec_Der(testCase,f_size,x0_size)
			% Robustesse face a l'orientation des vecteurs (rangee ou 
			% colonne) pour les arguments d'entrees (sauf pour la matrice jacobienne) 
			
			jac_fct	=	@(x) [5*0.1*x(2)*cos(0.1*x(1)*x(2)) , 5*0.1*x(1)*cos(0.1*x(1)*x(2)) , -1 ;...
							  2*x(1) , 2*x(2) , 2*x(3) ; ...
							  1 , -1 , -1];
					
			[app , err] = newton_ND_avec_der(f_size, jac_fct, x0_size , 20 , 1e-12);
			
			nb_iter		=	length(err);
			
			verifySize(testCase,app,[3,nb_iter]);
			verifySize(testCase,err,[1,nb_iter]);
			
		end
		
		function testBissecErrorFct(testCase)
			% Verification message erreur pour fct 
			verifyError(testCase,@() bissec(@fct_not_present, 1, 2, 200, 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() bissec(@fct_error, 1, 2, 200, 1e-12),"MATLAB:UndefinedFunction");
		end
		
		function testSecanteErrorFct(testCase)
			% Verification message erreur pour fct 
			verifyError(testCase,@() secante(@fct_not_present, 1, 2, 200, 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() secante(@fct_error, 1, 2, 200, 1e-12),"MATLAB:UndefinedFunction");
		end
		
		function testNewton_1DErrorFct(testCase)
			% Verification message erreur pour fct
			verifyError(testCase,@() newton_1D(@fct_not_present, @(x) 2*x, 1, 20, 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_1D(@(x) x.^2-10, @fct_not_present, 1, 20, 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_1D(@fct_error, @(x) 2*x, 1, 20, 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_1D(@(x) x.^2-10, @fct_error, 1, 20, 1e-12),"MATLAB:UndefinedFunction");
		end
		
		function testPts_FixesErrorFct(testCase)
			% Verification message erreur pour fct
			verifyError(testCase,@() pts_fixes(@fct_not_present, 1, 20, 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() pts_fixes(@fct_error, 1, 20, 1e-12),"MATLAB:UndefinedFunction");
		end
		
		function testNewton_ND_Avec_DerErrorFct(testCase)
			% Verification message erreur pour fct
			verifyError(testCase,@() newton_ND_avec_der(@f_not_present, @(x) [2*x(1),0;0,2*x(2)], [1;1] , 20 , 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_ND_avec_der(@(x) [x(1).^2-10;x(2).^2-11],@f_not_present, [1;1] , 20 , 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_ND_avec_der(@fct_error, @(x) [2*x(1),0;0,2*x(2)], [1;1] , 20 , 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_ND_avec_der(@(x) [x(1).^2-10;x(2).^2-11],@fct_error, [1;1] , 20 , 1e-12),"MATLAB:UndefinedFunction");
		end	
		
		function testNewton_ND_Sans_DerErrorFct(testCase)
			% Verification message erreur pour fct
			verifyError(testCase,@() newton_ND_sans_der(@f_not_present, 1e-3, [1;1] , 20 , 1e-12),"MATLAB:UndefinedFunction");
			verifyError(testCase,@() newton_ND_sans_der(@fct_error, 1e-3, [1;1] , 20 , 1e-12),"MATLAB:UndefinedFunction");
		end	
	end
	
end

function [ordre,ordre_app] = order_computation_nl(erreur,varargin)
% Approximation de l'ordre pour les methodes de resolution de problemes
% non-lineaires

	if nargin>2
		error("Il ne peut y avoir qu'un deuxieme argument, la tolerance specifiee.")
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
		warning("La zone asymptotique n'est pas tres grande")
	elseif any(gradient(ind_stable_region)~=1)
		warning("La zone asymptotique est brisee")
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
		error("Il ne peut y avoir qu'un troisieme argument, la tolerance specifiee.")
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
		warning("La zone asymptotique n'est pas tres grande")
	elseif any(gradient(ind_stable_region)~=1)
		warning("La zone asymptotique est brisee")
	end
	
	taux = mean(taux_app(ind_stable_region));
end

function [f] = fct_error(x)
	f = a*x;
end