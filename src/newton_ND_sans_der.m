function [approx , err_abs] = newton_ND_sans_der(F , x0 , nb_it_max , tol_rel)
% NEWTON_ND_SANS_DER	Méthode de Newton pour la résolution de F(x) = 0, pour F: R^n -> R^n
%
% Syntaxe: [approx , err_abs] = newton_ND_sans_der(F , x0 , nb_it_max , tol_rel)
%
% Argument d'entrée
%	F			-	String ou fonction handle spécifiant la fonction
%					non-linéaire (F: R^n -> R^n)
%	x0			-	Approximation initiale (x0 in R^n)
%	nb_it_max	-	Nombre maximum d'itérations 
%	tol			-	Tolérance sur l'approximation de l'erreur relative
%
% Arguments de sortie
%	approx		-	Vecteur colonne de taille nb_iter contenant les 
%					itérations
%	err_abs		-	Vecteur colonne de dimension nb_iter contenant les
%					erreurs absolues
%
% Exemples d'appel
%	[ approx , err_abs ] = newton_ND_sans_der( 'my_sys_nl' , [1,1] , 20 , 1e-9 )


%%  Vérification de la fonction contenant les dérivées
if isa(F,'char')
	fct		=	str2func(F);
elseif isa(F,'function_handle')
	fct		=	F;
else
	error('L''argument f n''est pas un string ni un function_handle')
end


%% Vérification du nb de composantes des conditions initiales et de f
if ~isnumeric(x0) || ~isvector(x0)
	error('L''approximation initiale x0 n''est pas un vecteur')
end

try 
	fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction F n''est pas dans le répertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction F ne retourne pas un vecteur')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction F doit prendre seulement 1 argument en entrée')
	else 
		rethrow(ME)
	end
end

taille	=	length(x0);

if ~isnumeric(fct(x0)) || ~isvector(fct(x0)) || (length(fct(x0))~=taille)
	error(['Le fonction F ne retourne pas un vecteur de même taille que x0. ',...
		'x0 est de taille %d alors que F(x0) est de taille %d.'],taille,length(fct(x0)))
end


%% Initialisation des matrices app et err
app			=	nan(nb_it_max,taille);
app(1,:)	=	x0;
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% Méthode de Newton
for t=1:nb_it_max-1
	
	app_jac		=	app_jacobienne(fct,app(t,:));
	delta_x		=	app_jac\-reshape(fct(app(t,:)),taille,1);
	app(t+1,:)	=	app(t,:) + delta_x';
	
	if min(abs(eig(app_jac))) == 0  
		warning(['La matrice jacobienne de f à l''itération %d est singulière 0.\n',...
					'Arrêt de l''algorithme'],t)
		break
	end
	
	err_rel(t)	=	norm(app(t+1,:)-app(t,:))/(norm(abs(app(t+1,:))) + eps);
	if (err_rel(t) <= tol_rel) || (norm(fct(app(t+1,:))) == 0)
		arret	=	true;
		break
	end
end

nb_it	=	t+1;
approx	=	app(1:nb_it,:);
err_abs	=	inf(nb_it,1);

if arret
	for t=1:nb_it
		err_abs(t)	=	norm(approx(end,:) - approx(t,:));
	end
else
	warning('La méthode de Newton n''a pas convergée')
end


end


function [app_finale] = app_jacobienne(f,x0)

	taille	=	length(x0);
	if min(x0) == 0
		h_init	=	1e-6;
	else
		h_init	=	1e-3 * min(x0);
	end
	h		=	h_init./(2.^(0:1));
	app		=	cell(2,1);
	
	for t=1:length(h)
		app{t}			=	zeros(taille,taille);
		for d=1:taille
			delta_h		=	zeros(size(x0));
			delta_h(d)	=	h(t);
			app{t}(:,d)	=	(f(x0+delta_h) - f(x0-delta_h))/(2*h(t));
		end
	end
	
	app_finale	=	(2^2*app{2} - app{1})/(2^2-1);
end