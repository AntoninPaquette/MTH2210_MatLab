function tests = interpolationTest
	tests = functiontests(localfunctions);
end

%% Vérification de la fonction lagrange

function testLagrangeSanity1(testCase)
	fct			=	@(x) cos(x) .* sin(x).^2;
	xi			=	[-10,-2,4,pi,-exp(1),0];
	yi			=	fct(xi);
	y_inter		=	lagrange(xi,yi,xi(5));
	y_exacte	=	fct(xi(5));
	err_rel		=	norm(y_exacte - y_inter)/norm(y_exacte);
	
	verifyLessThan(testCase,err_rel,1e-14);
end

function testLagrangeExactInter1(testCase)
% Interpolation exacte d'un polynome de degré 4

	fct			=	@(x) (x-1).*(x+2.5).*(x-pi).*(x+11).^2;
	xi			=	[-10,-2,4,pi,-exp(1),0];
	yi			=	fct(xi);
	xfin		=	linspace(min(xi),max(xi),1000);
	y_inter		=	lagrange(xi,yi,xfin);
	y_exacte	=	fct(xfin);
	err_rel		=	norm(y_exacte - y_inter)/norm(y_exacte);
	
	verifyLessThan(testCase,err_rel,1e-14);
end

function testLagrangeExactInter2(testCase)
% Interpolation exacte d'un polynome de degré 9

	fct			=	@(x) (x-2).^2 .* (x+5).^3 .* (x-exp(1)).^4;
	xi			=	linspace(-2,2,10);
	yi			=	fct(xi);
	xfin		=	linspace(min(xi),max(xi),1000);
	y_inter		=	lagrange(xi,yi,xfin);
	y_exacte	=	fct(xfin);
	err_rel		=	norm(y_exacte - y_inter)/norm(y_exacte);
	
	verifyLessThan(testCase,err_rel,1e-14);
end

function testLagrangeOrder(testCase)
% Ordre de convergence de l'intrerpolation de Lagrange avec les abscisses 
% de Chebyshev

	fct			=	@(x) x.^3 .* exp(2*x);
	degre		=	3;
	nb_pts		=	degre + 1;
	nb_loop		=	10;
	x_interet	=	(1/3)^nb_loop;
	erreur		=	nan(nb_loop,1);
	erreur_rel	=	nan(nb_loop,1);
	y_exact		=	fct(x_interet);
	for t=1:nb_loop
		a			=	0;
		b			=	(1/2)^(t-1);
		k			=	1:nb_pts;
		xi			=	1/2*(a+b) + 1/2*(b-a)*cos((2*k-1)*pi/(2*nb_pts));
		yi			=	fct(xi);
		y_inter		=	lagrange(xi,yi,x_interet);
		erreur(t)	=	abs(y_exact - y_inter);
	end
	
	tol = 0.2;
	[ordre,~] = order_computation(erreur,2,tol);
	
	verifyLessThan(testCase,abs(ordre-(degre+1)),tol);
end

%% Vérification de la fonction splinec

function testSplineExact1(testCase)
% Interpolation exacte avec une spline pour b.c. naturelle/courbure
% prescrite

	fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
% 	d_fct	=	@(x) -12*x.^2 + 2*pi*x + 11;
	d2_fct	=	@(x) -24*x + 2*pi;
	a		=	pi/12;
	b		=	4;
	nb_pts	=	10;
	xi	=	linspace(a,b,nb_pts);
	yi	=	fct(xi);
	x	=	linspace(a,b,1000);
	y_exacte	=	fct(x);

	[ Sx ] = splinec( xi , yi , x , [1,2] , [0,d2_fct(b)]);

	err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
	verifyLessThan(testCase,err_rel,1e-14);
end

function testSplineExact2(testCase)
% Interpolation exacte avec une spline pour b.c. courbure
% prescrite/naturelle

	fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
% 	d_fct	=	@(x) -12*x.^2 + 2*pi*x + 11;
	d2_fct	=	@(x) -24*x + 2*pi;
	a		=	-10;
	b		=	pi/12;
	nb_pts	=	10;
	xi	=	linspace(a,b,nb_pts);
	yi	=	fct(xi);
	x	=	linspace(a,b,1000);
	y_exacte	=	fct(x);

	[ Sx ] = splinec( xi , yi , x , [2,1] , [d2_fct(a),0]);

	err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
	verifyLessThan(testCase,err_rel,1e-14);
end

function testSplineExact3(testCase)
% Interpolation exacte avec une spline pour b.c. courbure prescrite

	fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
% 	d_fct	=	@(x) -12*x.^2 + 2*pi*x + 11;
	d2_fct	=	@(x) -24*x + 2*pi;
	a		=	-5;
	b		=	4;
	nb_pts	=	10;
	xi	=	linspace(a,b,nb_pts);
	yi	=	fct(xi);
	x	=	linspace(a,b,1000);
	y_exacte	=	fct(x);

	[ Sx ] = splinec( xi , yi , x , [2,2] , [d2_fct(a),d2_fct(b)]);

	err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
	verifyLessThan(testCase,err_rel,1e-14);
end

function testSplineExact4(testCase)
% Interpolation exacte avec une spline pour b.c. pente prescrite

	fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
	d_fct	=	@(x) -12*x.^2 + 2*pi*x + 11;
% 	d2_fct	=	@(x) -24*x + 2*pi;
	a		=	-5;
	b		=	4;
	nb_pts	=	3;
	xi	=	linspace(a,b,nb_pts);
	yi	=	fct(xi);
	x	=	linspace(a,b,1000);
	y_exacte	=	fct(x);

	[ Sx ] = splinec( xi , yi , x , [4,4] , [d_fct(a),d_fct(b)]);

	err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
	verifyLessThan(testCase,err_rel,1e-14);
end

function testSplineExact5(testCase)
% Interpolation exacte avec une spline pour abscisse non uniformément
% distribuées

	fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
	d_fct	=	@(x) -12*x.^2 + 2*pi*x + 11;
% 	d2_fct	=	@(x) -24*x + 2*pi;
	a		=	-5;
	b		=	4;
	nb_pts	=	3;
	xi	=	[linspace(a,b,nb_pts),linspace(b+1,b+2,nb_pts)];
	yi	=	fct(xi);
	x	=	linspace(a,b+2,1000);
	y_exacte	=	fct(x);

	[ Sx ] = splinec( xi , yi , x , [4,4] , [d_fct(a),d_fct(b+2)]);

	err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
	verifyLessThan(testCase,err_rel,1e-14);
end

function testSplineConstantCourbure(testCase)
% Interpolation exacte avec une spline pour b.c. courbure constante

	xi	=	[0,1,3,4];
	yi	=	spline_courbure_constante(xi);
	x	=	linspace(1,4,1000);
	y_exacte	=	spline_courbure_constante(x);

	[ Sx ] = splinec( xi , yi , x , [3,3] , [NaN,NaN]);

	err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
	verifyLessThan(testCase,err_rel,1e-14);
end

function [px] = spline_courbure_constante(x)
	if min(x)<0 || max(x)>4
		error("Pas dans le domaine de la fonction")
	end
	px =	(x>=0 & x<1) .* (x.^2) + (x>=1 & x<3) .* (-x.^3 + 4*x.^2 - 3*x + 1) + ...
			(x>=3 & x<=4) .* (-5*x.^2 + 24*x - 26);
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
	stable_region		=	(ordre_app>0) & (abs(gradient(ordre_app))<tol);
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