function [approx , err_abs] = secante(f , x0 , x1 , nb_it_max , tol_rel)
% SECANTE	Méthode de la sécante pour la résolution f(x) = 0
%			pour f: R -> R
%
% Syntaxe: [approx , err_abs] = secante(f , x0 , x1 , nb_it_max , tol_rel)
%
% Argument d'entrée
%	f			-	String ou fonction handle spécifiant la fonction
%					non-linéaire
%	x0			-	1ère approximation initiale 
%	x1			-	2ème approximation initiale 
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
%	[ approx , err_abs ] = secante( @(x) x.^2-10 , 3 , 4 , 20 , 1e-9 )


%%  Vérification de la fonction contenant les dérivées
if isa(f,'char')
	fct		=	str2func(f);
elseif isa(f,'function_handle')
	fct		=	f;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

%% Vérification du nb de composantes des conditions initiales et de f
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
elseif ~isnumeric(x1) || ~isscalar(x1)
	error('L''approximation initiale x1 n''est pas un scalaire')
end

try 
	fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction f n''est pas dans le répertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction f ne retourne pas un scalaire')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction f doit prendre seulement 1 argument en entrée')
	else 
		rethrow(ME)
	end
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0)) 
	error('Le fonction f ne retourne pas un scalaire')
end


%% Initialisation des matrices app et err
app			=	nan(nb_it_max,1);
app(1)		=	x0;
app(2)		=	x1;
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% Méthode de la sécante
for t=2:nb_it_max-1
	
	app(t+1)	=	app(t) - fct(app(t)) * (app(t) - app(t-1))/...
										   (fct(app(t)) - fct(app(t-1)));
									   
	if abs(fct(app(t)) - fct(app(t-1))) == 0  
		warning(['L''approximation de la dérivée de f à l''aide des ',...
			   'points x_i=%6.5e et x_{i-1}=%6.5e ',...
			   'est exactement 0.\nArrêt de l''algorithme\n'],app(t),app(t-1))
		break
	end
	
	err_rel(t)	=	abs(app(t+1)-app(t))/(abs(app(t+1)) + eps);
	if (err_rel(t) <= tol_rel) || (fct(app(t+1)) == 0)
		arret	=	true;
		break
	end
end

nb_it	=	t+1;
approx	=	app(1:nb_it);
err_abs		=	inf(nb_it,1);

if arret
	err_abs		=	abs(approx(end) - approx);
else
	warning('La méthode de la sécante n''a pas convergée')
end


end