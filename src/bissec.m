function [approx , err_abs] = bissec(f , x0 , x1 , nb_it_max , tol_rel)
% BISSEC	Méthode de la bissection pour la résolution f(x) = 0
%			pour f: R -> R
%
% Syntaxe: [approx , err_abs] = bissec(f , x0 , x1 , nb_it_max , tol_rel)
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
%	[ approx , err_abs ] = bissec( @(x) x.^2-10 , 3 , 4 , 100 , 1e-9 )


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
elseif fct(x0)*fct(x1)>0
	warning('La condition f(x0)*f(x1)<0 n''est pas respectée.\nArrêt de l''algorithme%\n',[])
	approx	=	[x0;x1];
	err_abs =	inf(2,1);
	return
elseif fct(x0)==0
	approx	=	x0;
	err_abs	=	0;
	return
elseif fct(x1)==0
	approx	=	x1;
	err_abs	=	0;
	return
end

%% Initialisation des matrices app et err
app			=	nan(nb_it_max,1);
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% Méthode de la bissection
for t=1:nb_it_max-1
	
	if t==1
		x_gauche	=	min([x0,x1]);
		x_droite	=	max([x0,x1]);
	else
		if f_gauche*f_milieu < 0
			x_droite	=	x_milieu;
		elseif f_droite*f_milieu < 0
			x_gauche	=	x_milieu;
		else
			warning('Problème avec la fonction f.\nArrêt de l''algorithme%\n',[])
			break
		end
	end
	
	x_milieu	=	(x_gauche + x_droite)/2;
	app(t)		=	x_milieu;
	
	if t==1
		if fct(app(t)) == 0
			arret	=	true;
			break
		end
	else
		err_rel(t-1)	=	abs(app(t)-app(t-1))/(abs(app(t)) + eps);
		
		if (err_rel(t-1) <= tol_rel) || (fct(app(t)) == 0)
			arret	=	true;
			break
		end
	end

	f_gauche	=	fct(x_gauche);
	f_droite	=	fct(x_droite);
	f_milieu	=	fct(x_milieu);
	
end

nb_it	=	t;
approx	=	app(1:nb_it);
err_abs	=	inf(nb_it,1);

if arret
	err_abs		=	abs(approx(end) - approx);
else
	warning('La méthode de la bissection n''a pas convergée')
end


end