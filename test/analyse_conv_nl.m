%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script effectuant l'analyse de convergence des méthodes de résolution de
% problèmes non-linéaires
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

clc
close all

alpha	=	(1+sqrt(5))/2;

%% Vérification des ordre de convergence pour bissect, secante, newton1D
% Premier problème
fct1	=	@(x) x^2 - 10;
d_fct1	=	@(x) 2*x;

x0	=	2;
x1	=	5;


[approx_bis1 , err_bis1] = bissec(fct1 , x0 , x1 , 200 , 1e-12);
[approx_sec1 , err_sec1] = secante(fct1 , x0 , x1 , 200 , 1e-12);
[approx_new1 , err_new1] = newton_1D(fct1 , d_fct1, x0 , 200 , 1e-12);

err_ex_bis1		=	abs(approx_bis1-sqrt(10));
err_ex_sec1		=	abs(approx_sec1-sqrt(10));
err_ex_new1		=	abs(approx_new1-sqrt(10));

order_computation_bissect(err_sec1)
order_computation_nl(err_sec1,0.2)

% Calcul des ratios

ratio1_sec_fct1			=	err_ex_sec1(2:end)./err_ex_sec1(1:end-1);
ratio_alpha_sec_fct1	=	err_ex_sec1(2:end)./err_ex_sec1(1:end-1).^alpha;
ratio2_sec_fct1			=	err_ex_sec1(2:end)./err_ex_sec1(1:end-1).^2;

ratio1_new_fct1			=	err_ex_new1(2:end)./err_ex_new1(1:end-1);
ratio_alpha_new_fct1	=	err_ex_new1(2:end)./err_ex_new1(1:end-1).^alpha;
ratio2_new_fct1			=	err_ex_new1(2:end)./err_ex_new1(1:end-1).^2;

figure
loglog(err_ex_bis1(1:end-1),err_ex_bis1(2:end))
hold on
loglog(err_ex_sec1(1:end-1),err_ex_sec1(2:end))
loglog(err_ex_new1(1:end-1),err_ex_new1(2:end))
legend('Bissection','Secante','Newton')

figure
semilogy(1:length(err_ex_bis1),err_ex_bis1)
hold on
semilogy(1:length(err_ex_sec1),err_ex_sec1)
semilogy(1:length(err_ex_new1),err_ex_new1)
legend('Bissection','Secante','Newton')
%% Vérification des ordre de convergence pour bissect, secante, newton1D
% Deuxième problème

fct2	=	@(x) exp(x) - x^3;
d_fct2	=	@(x) exp(x) - 3*x^2;

x0	=	1.5;
x1	=	2.5;


[approx_bis2 , err_bis2] = bissec(fct2 , x0 , x1 , 200 , 1e-12);
[approx_sec2 , err_sec2] = secante(fct2 , x0 , x1 , 200 , 1e-12);
[approx_new2 , err_new2] = newton_1D(fct2 , d_fct2, x0 , 200 , 1e-12);

% Calcul des ratios

ratio1_sec_fct2			=	err_sec2(2:end)./err_sec2(1:end-1);
ratio_alpha_sec_fct2	=	err_sec2(2:end)./err_sec2(1:end-1).^alpha;
ratio2_sec_fct2			=	err_sec2(2:end)./err_sec2(1:end-1).^2;

ratio1_new_fct2			=	err_new2(2:end)./err_new2(1:end-1);
ratio_alpha_new_fct2	=	err_new2(2:end)./err_new2(1:end-1).^alpha;
ratio2_new_fct2			=	err_new2(2:end)./err_new2(1:end-1).^2;

figure
semilogy(1:length(err_bis2),err_bis2)
hold on
semilogy(1:length(err_sec2),err_sec2)
semilogy(1:length(err_new2),err_new2)
legend('Bissection','Secante','Newton')


% Vérification de Newton pour racines multiples

fct3	=	@(x) x*sin(x)^2;
d_fct3	=	@(x) sin(x)^2 + 2*x*sin(x)*cos(x);

[approx_new3 , err_new3] = newton_1D(fct3 , d_fct3, 1 , 200 , 1e-12);
[approx_new4 , err_new4] = newton_1D(fct3 , d_fct3, 3 , 200 , 1e-12);

err_ex_new3		=	abs(approx_new3-0);
err_ex_new4		=	abs(approx_new4-pi);

ratio1_new_fct3		=	err_ex_new3(2:end)./err_ex_new3(1:end-1); % Racine de multiplicité 3
ratio1_new_fct4		=	err_ex_new4(2:end)./err_ex_new4(1:end-1); % Racine de multiplicité 2


%% Vérification de la méthode des points-fixes pour ordre 1

fct4	=	@(x) -x^2/10 + x + 1;
taux	=	-2*sqrt(10)/10 + 1;

[approx_ptf5 , err_ptf5] = pts_fixes(fct4 , 1 , 200 , 1e-12);

err_ex_ptf5		=	abs(approx_ptf5 - sqrt(10));

ratio1_ptf5		=	err_ex_ptf5(2:end)./err_ex_ptf5(1:end-1);


%% Vérification de la méthode des points-fixes pour ordre 2

fct5	=	@(x) -x^2/6 + x + 9/6;
taux2	=	-1/6;

[approx_ptf6 , err_ptf6] = pts_fixes(fct5 , 1 , 200 , 1e-12);

err_ex_ptf6		=	abs(approx_ptf6 - 3);

ratio1_ptf6		=	err_ex_ptf6(2:end)./err_ex_ptf6(1:end-1);
ratio2_ptf6		=	err_ex_ptf6(2:end)./err_ex_ptf6(1:end-1).^2;



%% Vérification des fonctions newton_ND_avec_der et newton_ND_sans_der

fct6	=	@(x) [5*sin(0.1*x(1)*x(2)) - x(3) , ...
				  x(1)^2 + x(2)^2 + x(3)^2 - 9 , ...
				  x(1) - x(2) - x(3) - 1];
			  
jac_fct6	=	@(x) [5*0.1*x(2)*cos(0.1*x(1)*x(2)) , 5*0.1*x(1)*cos(0.1*x(1)*x(2)) , -1 ;...
				 2*x(1) , 2*x(2) , 2*x(3) ; ...
				 1 , -1 , -1];

x0 = [1 , 1 , 1];
			 
[approx_sans_6 , err_sans_6] = newton_ND_sans_der(fct6 , x0 , 100 , 1e-12);
[approx_avec_6 , err_avec_6] = newton_ND_avec_der(fct6 , jac_fct6 , x0 , 100 , 1e-12);

ratio1_sans_6	=	err_sans_6(2:end)./err_sans_6(1:end-1);
ratio_alpha_sans_6	=	err_sans_6(2:end)./err_sans_6(1:end-1).^alpha;
ratio2_sans_6	=	err_sans_6(2:end)./err_sans_6(1:end-1).^2;

ratio1_avec_6	=	err_avec_6(2:end)./err_avec_6(1:end-1);
ratio_alpha_avec_6	=	err_avec_6(2:end)./err_avec_6(1:end-1).^alpha;
ratio2_avec_6	=	err_avec_6(2:end)./err_avec_6(1:end-1).^2;


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
		error("Il ne pas y avoir de zone asymptotique")
	elseif length(ind_stable_region) < 2
		warning("La zone asymptotique n'est pas très grande")
	elseif any(gradient(ind_stable_region)~=1)
		warning("La zone asymptotique est brisée")
	end
	
	ordre = mean(ordre_app(ind_stable_region));
end

function [ordre] = order_computation_bissect(erreur)
	
	% Least-square fit
	nb_iter	=	length(erreur);
	A		=	[ones(nb_iter-2,1) reshape(log(erreur(1:end-2)),[],1)];
	b		=	reshape(log(erreur(2:end-1)),[],1);
	coeff	=	A\b;
	
	ordre	=	coeff(2);

end