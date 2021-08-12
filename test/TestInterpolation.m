classdef TestInterpolation < matlab.unittest.TestCase
	
	properties (TestParameter)
		% Interpolation exacte de polynomes
		fct_degre	=	{@(x) (x-1).*(x+2.5).*(x-pi).*(x+11), @(x) (x-2).^2 .* (x+5).^3 .* (x-exp(1)).^4} 
		degre		=	{4,9}
		
		degre_for_order =	{1,2,3,4,5}
				
		% Spline décrite par la fonction spline_example ci-dessous avec
		% bc1 en x=0 et bc2 en x=4. Les conditions frontière sont:
		% S'(0)=0 , S''(0) = 2 , S'(4) = -16 et S''(4) = -10.
		spline_type_bc1	=	{2,2,2,3,3,3,4,4,4}
		spline_type_bc2	=	{2,3,4,2,3,4,2,3,4}
		spline_bc1		=	{2,2,2,NaN,NaN,NaN,0,0,0}
		spline_bc2		=	{-10,NaN,-16,-10,NaN,-16,-10,NaN,-16}
			
	end
	
	methods (Test, ParameterCombination = 'sequential')
		function testLagrangeExact(testCase,fct_degre,degre)
			% Interpolation exacte de polynome
			
			xi	=	linspace(-5,5,degre+1);
			yi	=	fct_degre(xi);
			x_fin	=	linspace(min(xi),max(xi),1000);
			y_inter	=	lagrange(xi,yi,x_fin);
			y_exacte	=	fct_degre(x_fin);
			err_rel		=	norm(y_exacte - y_inter)/norm(y_exacte);
			
			verifyLessThan(testCase,err_rel,1e-14);
		end
		
		function testSplineExact1(testCase,spline_type_bc1,spline_type_bc2,spline_bc1,spline_bc2)
			% Interpolation exacte de la spline décrite dans la fonction 
			% spline_example avec toute les combinaisons de condition
			% frontière courbure prescrite (2), courbure constante (3) et
			% pente prescrite (4)
			
			xi	=	[0,1,3,4];
			yi	=	spline_example2(xi);
			x	=	linspace(1,4,1000);
			y_exacte	=	spline_example(x);
			
			[ Sx ] = splinec( xi , yi , x , [spline_type_bc1,spline_type_bc2] , [spline_bc1,spline_bc2]);
			
			err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
			verifyLessThan(testCase,err_rel,1e-14);
		end
	end
	
	methods (Test)
		function testLagrangeOrder(testCase, degre_for_order)
			% Ordre de convergence de l'intrerpolation de Lagrange avec 
			% les abscisses de Chebyshev
			
			fct			=	@(x) exp(2*x);
			nb_pts		=	degre_for_order + 1;
			nb_loop		=	10;
			x_interet	=	(1/3)^nb_loop;
			erreur		=	nan(nb_loop,1);
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
			
			[ordre,~] = order_computation(erreur,2,0.2);
			
			verifyLessThan(testCase,abs(ordre-(degre_for_order+1)),0.2);
		end
		
		function testSplineNatCourbure(testCase)
			% Spline avec condition naturelle et courbure prescrite
			
			fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
			d2_fct	=	@(x) -24*x + 2*pi;
			a		=	pi/12;
			b		=	4;
			nb_pts	=	10;
			xi	=	linspace(a,b,nb_pts);
			yi	=	fct(xi);
			x	=	linspace(a,b,1000);
			y_exacte	=	fct(x);
			
			[ Sx ] = splinec( xi , yi , x , [1,2] , [NaN,d2_fct(b)]);
			
			err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
			verifyLessThan(testCase,err_rel,1e-14);
		end
		
		function SplinePenteNat(testCase)
			% Spline avec pente prescrite et condition naturelle

			fct		=	@(x) -4*x.^3 + pi*x.^2 + 11*x - exp(1);
			d_fct	=	@(x) -12*x.^2 + 2*pi*x + 11;
			a		=	-10;
			b		=	pi/12;
			nb_pts	=	10;
			xi	=	linspace(a,b,nb_pts);
			yi	=	fct(xi);
			x	=	linspace(a,b,1000);
			y_exacte	=	fct(x);
			
			[ Sx ] = splinec( xi , yi , x , [4,1] , [d_fct(a),NaN]);
			
			err_rel	=	norm(y_exacte - Sx)/norm(y_exacte);
			verifyLessThan(testCase,err_rel,1e-14);
		end
	end
end

function [px] = spline_example(x)
	if min(x)<0 || max(x)>4
		error("Pas dans le domaine de la fonction")
	end
	px =	(x>=0 & x<1) .* (x.^2) + (x>=1 & x<3) .* (-x.^3 + 4*x.^2 - 3*x + 1) + ...
			(x>=3 & x<=4) .* (-5*x.^2 + 24*x - 26);
end

function [px] = spline_example2(x)
	if min(x)<0 || max(x)>4
		error("Pas dans le domaine de la fonction")
	end
	px =	(x>=0 & x<1) .* (x.^2) + (x>=1 & x<3) .* (-x.^3 + 4*x.^2 - 3*x + 1) + ...
			(x>=3 & x<=4) .* (-5*x.^2 + 24*x - 26);
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
