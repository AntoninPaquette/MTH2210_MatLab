function tests = nonLineaireTest()
	tests = functiontests(localfunctions);
end


function setup(testCase)
	testCase.TestData.alpha		=	(1+sqrt(5))/2;
	testCase.TestData.tol_alg	=	1e-12;
	
	testCase.TestData.f1		=	@(x) x^2 - 10;
	testCase.TestData.d_fct1	=	@(x) 2*x;
	testCase.TestData.r1		=	sqrt(10);
	testCase.TestData.x0_f1		=	2;
	testCase.TestData.x1_f1		=	5;
	
	testCase.TestData.f2		=	@(x) exp(x) - x^3;
	testCase.TestData.d_fct2	=	@(x) exp(x) - 3*x^2;
	testCase.TestData.x0_f2		=	1.5;
	testCase.TestData.x1_f2		=	2.5;

end

function bissectOrder1Test(testCase)

	[approx, err] = bissec(testCase.TestData.f1, testCase.TestData.x0_f1, ...
									testCase.TestData.x1_f1 , 200 , testCase.TestData.tol_alg);
	
	[ordre] = order_computation_bissect(err);							
	
	verifyLessThan(testCase,abs(approx(end)-testCase.TestData.r1),testCase.TestData.tol_alg);
	verifyLessThan(testCase,abs(ordre-1),0.2);
end


function [ordre,ordre_app] = order_computation_nl(erreur,varargin)

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
	
	% Least-square fit
	nb_iter	=	length(erreur);
	A		=	[ones(nb_iter-2,1) reshape(log(erreur(1:end-2)),[],1)];
	b		=	reshape(log(erreur(2:end-1)),[],1);
	coeff	=	A\b;
	
	ordre	=	coeff(2);

end

