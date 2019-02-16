function [approx , err_abs] = newton_1D(f , df , x0 , nb_it_max , tol_rel , file_name)
% NEWTON_1D	Méthode de Newton pour la résolution f(x) = 0,
%			pour f: R -> R
%
% Syntaxe: [approx , err_abs] = newton_1D(f , df , x0 , nb_it_max ,tol_rel)
%
% Argument d'entrée
%	f			-	String ou fonction handle spécifiant la fonction
%					non-linéaire
%	df			-	String ou function handle spécifiant la dérivée de f
%	x0			-	Approximation initiale
%	nb_it_max	-	Nombre maximum d'itï¿½rations
%	tol_rel		-	Tolérance sur l'approximation de l'erreur relative
%	file_name	-	(Optionnel) Nom du fichier (avec l'extension .txt) dans
%					lequel sera	écrit les résultats de l'algorithme
%
% Arguments de sortie
%	approx		-	Vecteur colonne de taille nb_iter contenant les
%					itérations
%	err_abs		-	Vecteur colonne de dimension nb_iter contenant les
%					erreurs absolues
%
% Exemples d'appel
%	[ approx , err_abs ] = newton_1D( 'my_fct_nl' , 'my_dfct_nl' , 3 , 20 , 1e-9 , 'resul_newton.txt')
%	[ approx , err_abs ] = newton_1D( @(x) x.^2-10 , @(x) 2*x , 3 , 20 , 1e-9 , , 'resul_newton.txt')



%%  Vérification de la fonction contenant les dérivées
if isa(f,'char')
	fct			=	str2func(f);
	is_fct_file =	true;
elseif isa(f,'function_handle')
	fct			=	f;
	is_fct_file =	false;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

if isa(df,'char')
	d_fct			=	str2func(df);
	is_dfct_file	=	true;
elseif isa(df,'function_handle')
	d_fct			=	df;
	is_dfct_file	=	false;
else
	error('L''argument df n''est pas un string ni un function_handle')
end

%% Vérification du nb de composantes des conditions initiales et de f
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
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

try
	d_fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction df n''est pas dans le répertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction df ne retourne pas un scalaire')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction df doit prendre seulement 1 argument en entrée')
	else
		rethrow(ME)
	end
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0))
	error('Le fonction f ne retourne pas un scalaire')
elseif ~isnumeric(d_fct(x0)) || ~isscalar(d_fct(x0))
	error('Le fonction df ne retourne pas un scalaire')
elseif ~check_derivative(fct,d_fct,x0)
	warning('Il semble y avoir une erreur avec la dérivées')
end

%% Initialisation des matrices app et err
app			=	nan(nb_it_max,1);
app(1)		=	x0;
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% Méthode de Newton
for t=1:nb_it_max-1

	app(t+1)	=	app(t) - fct(app(t))/d_fct(app(t));

	if abs(d_fct(app(t))) == 0
		warning(['La dérivée de f au point x=%6.5e est exactement 0.\n',...
					'Arrêt de l''algorithme'],app(t))
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
	warning('La méthode de Newton n''a pas convergée')
end

% Écriture des résultats si fichier passé en argument
if nargin == 6
	output_results(file_name , fct , is_fct_file , d_fct, is_dfct_file, ...
				nb_it_max , tol_rel , x0 , approx , err_abs , arret)
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

function [] = output_results(file_name , fct , is_fct_file , d_fct, ...
			is_dfct_file,it_max , tol_rel , x0 , x , err , status)
						 
	fid		=	fopen(file_name,'w');
	fprintf(fid,'Algorithme de Newton avec dérivée\n\n');
	fprintf(fid,'Fonction dont on cherche les racines:\n');
	if is_fct_file
		fprintf(fid,'%s\n\n',fileread([func2str(fct),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(fct));
	end
	fprintf(fid,'Dérivée de la fonction dont on cherche les racines:\n');
	if is_dfct_file
		fprintf(fid,'%s\n\n',fileread([func2str(d_fct),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(d_fct));
	end
	fprintf(fid,'Arguments d''entrée:\n');
	fprintf(fid,'    - Nombre maximum d''itérations: %d\n',it_max);
	fprintf(fid,'    - Tolérance relative: %6.5e\n',tol_rel);
	fprintf(fid,'    - Approximation initiale x0: %16.15e\n',x0);
	
	if status
		fprintf(fid,'\nStatut: L''algorithme de Newton a convergé en %d itérations\n\n',length(x)-1);
	else
		fprintf(fid,'\nStatut: L''algorithme de Newton n''a pas convergé\n\n');
	end
	fprintf(fid,'#It            x              Erreur absolue\n');
	fprintf(fid,'%3d   %16.15e   %6.5e\n',[reshape(0:length(x)-1,1,[]);reshape(x,1,[]);reshape(err,1,[])]);
	fclose(fid);

end
