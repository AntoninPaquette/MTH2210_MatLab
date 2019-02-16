function [approx , err_abs] = pts_fixes(g , x0 , nb_it_max , tol_rel , file_name)
% PTS_FIXES	Méthode des points-fixes pour la résolution de g(x) = x
%			pour g: R -> R
%
% Syntaxe: [approx , err_abs] = pts_fixes(g , x0 , nb_it_max , tol_rel)
%
% Argument d'entrée
%	g			-	String ou fonction handle spécifiant la fonction
%					non-linéaire
%	x0			-	Approximation initiale 
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
%	[approx , err_abs] = pts_fixes(@(x) -x.^2/10+x+1 , 4 , 50 , 1e-12 , 'resul_pts_fixes.txt')



%%  Vérification de la fonction contenant les dérivées
if isa(g,'char')
	fct			=	str2func(g);
	is_fct_file =	true;
elseif isa(g,'function_handle')
	fct			=	g;
	is_fct_file =	false;
else
	error('L''argument g n''est pas un string ni un function_handle')
end

%% Vérification du nb de composantes des conditions initiales et de f
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
end

try 
	fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction g n''est pas dans le répertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction g ne retourne pas un scalaire')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction g doit prendre seulement 1 argument en entrée')
	else 
		rethrow(ME)
	end
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0)) 
	error('Le fonction g ne retourne pas un scalaire')
end

%% Initialisation des matrices app et err
app			=	nan(nb_it_max,1);
app(1)		=	x0;
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% Méthode des points fixes
for t=1:nb_it_max-1	
	app(t+1)	=	fct(app(t));
	err_rel(t)	=	abs(app(t+1)-app(t))/(abs(app(t+1)) + eps);
	if err_rel(t) <= tol_rel || fct(app(t+1)) == app(t+1)
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
	warning('La méthode des points-fixes n''a pas convergée')
end

if nargin == 5
	output_results(file_name , fct , is_fct_file , ...
				nb_it_max , tol_rel , x0 , approx , err_abs , arret)
end

end


function [] = output_results(file_name , fct , is_fct_file , ...
			it_max , tol_rel , x0 , x , err , status)
						 
	fid		=	fopen(file_name,'w');
	fprintf(fid,'Algorithme des points-fixes\n\n');
	fprintf(fid,'Fonction dont on cherche les points-fixes:\n');
	if is_fct_file
		fprintf(fid,'%s\n\n',fileread([func2str(fct),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(fct));
	end

	fprintf(fid,'Arguments d''entrée:\n');
	fprintf(fid,'    - Nombre maximum d''itérations: %d\n',it_max);
	fprintf(fid,'    - Tolérance relative: %6.5e\n',tol_rel);
	fprintf(fid,'    - Approximation initiale x0: %16.15e\n',x0);
	
	if status
		fprintf(fid,'\nStatut: L''algorithme des points-fixes a convergé en %d itérations\n\n',length(x)-1);
	else
		fprintf(fid,'\nStatut: L''algorithme des points-fixes n''a pas convergé\n\n');
	end
	fprintf(fid,'#It            x              Erreur absolue\n');
	fprintf(fid,'%3d   %16.15e   %6.5e\n',[reshape(0:length(x)-1,1,[]);reshape(x,1,[]);reshape(err,1,[])]);
	fclose(fid);

end