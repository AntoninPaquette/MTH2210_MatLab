function [approx , err_abs] = newton_1D(f , df , x0 , nb_it_max , tol_rel)
% NEWTON_1D	M�thode de Newton pour la r�solution f(x) = 0,
%			pour f: R -> R
%
% Syntaxe: [approx , err_abs] = newton_1D(f , df , x0 , nb_it_max ,tol_rel)
%
% Argument d'entr�e
%	f			-	String ou fonction handle sp�cifiant la fonction
%					non-lin�aire
%	df			-	String ou function handle sp�cifiant la d�riv�e de f
%	x0			-	Approximation initiale
%	nb_it_max	-	Nombre maximum d'it�rations
%	tol_rel			-	Tol�rance sur l'approximation de l'erreur relative
%
% Arguments de sortie
%	approx		-	Vecteur colonne de taille nb_iter contenant les
%					it�rations
%	err_abs		-	Vecteur colonne de dimension nb_iter contenant les
%					erreurs absolues
%
% Exemples d'appel
%	[ approx , err_abs ] = newton_1D( @(x) x.^2-10 , @(x) 2*x , 3 , 20 , 1e-9 )


%%  V�rification de la fonction contenant les d�riv�es
if isa(f,'char')
	fct		=	str2func(f);
elseif isa(f,'function_handle')
	fct		=	f;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

if isa(df,'char')
	d_fct		=	str2func(df);
elseif isa(df,'function_handle')
	d_fct		=	df;
else
	error('L''argument df n''est pas un string ni un function_handle')
end

%% V�rification du nb de composantes des conditions initiales et de f
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
end

try
	fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction f n''est pas dans le r�pertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction f ne retourne pas un scalaire')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction f doit prendre seulement 1 argument en entr�e')
	else
		rethrow(ME)
	end
end

try
	d_fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction df n''est pas dans le r�pertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction df ne retourne pas un scalaire')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction df doit prendre seulement 1 argument en entr�e')
	else
		rethrow(ME)
	end
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0))
	error('Le fonction f ne retourne pas un scalaire')
elseif ~isnumeric(d_fct(x0)) || ~isscalar(d_fct(x0))
	error('Le fonction df ne retourne pas un scalaire')
elseif ~check_derivative(fct,d_fct,x0)
	warning('Il semble y avoir une erreur avec la d�riv�es')
end

%% Initialisation des matrices app et err
app			=	nan(nb_it_max,1);
app(1)		=	x0;
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% M�thode de Newton
for t=1:nb_it_max-1

	app(t+1)	=	app(t) - fct(app(t))/d_fct(app(t));

	if abs(d_fct(app(t))) == 0
		warning(['La d�riv�e de f au point x=%6.5e est exactement 0.\n',...
					'Arr�t de l''algorithme'],app(t))
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
	warning('La m�thode de Newton n''a pas converg�e')
end


end


function [test] = check_derivative(f,df,x0)

	if x0 == 0
		h_init	=	1e-6;
	else
		h_init	=	1e-3 * abs(x0);
	end

	h		=	h_init./(2.^(0:4));
	erreur	=	nan(size(h));

	for t=1:length(h)
		app			=	(f(x0+h(t)) - f(x0-h(t)))/(2*h(t));
		erreur(t)	=	abs(df(x0) - app);
	end

	ordre	=	log(erreur(2:end) ./ erreur(1:end-1)) ./ ...
				log(h(2:end)      ./ h(1:end-1));

	test	=	(mean( abs(2 - ordre) ) <= 0.25) || (mean(erreur/(abs(df(x0))+eps)) <=1e-9);

end
