function [approx , err_abs] = secante(f , x0 , x1 , nb_it_max , tol_rel , file_name)
% SECANTE	Méthode de la sécante pour la résolution f(x) = 0
%			pour f: R -> R
%
% Syntaxe: [approx , err_abs] = secante(f , x0 , x1 , nb_it_max , tol_rel)
%
% Arguments d'entrée
%	f			-	String ou fonction handle spécifiant la fonction
%					non-linéaire
%	x0			-	1ère approximation initiale
%	x1			-	2ème approximation initiale
%	nb_it_max	-	Nombre maximum d'itérations
%	tol_rel		-	Tolérance sur l'approximation de l'erreur relative
%	file_name	-	(Optionnel) Nom du fichier (avec l'extension .txt) dans
%					lequel sera	écrit les résultats de l'algorithme
%
% Arguments de sortie
%	approx		-	Vecteur rangee de taille nb_iter contenant les
%					itérations
%	err_abs		-	Vecteur rangee de dimension nb_iter contenant les
%					erreurs absolues
%
% Exemples d'appel
%	[ approx , err_abs ] = secante( 'my_fct_nl' , 3 , 4 , 100 , 1e-9 , 'resul_secante.txt')
%	[ approx , err_abs ] = secante( @(x) x.^2-10 , 3 , 4 , 100 , 1e-9 , 'resul_secante.txt')


%%  Vérification de la fonction f
if isa(f,'char')
	fct		=	str2func(f);
	is_fct_file = true;
elseif isa(f,'function_handle')
	fct		=	f;
	is_fct_file = false;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

%% Vérification du nb de composantes de l'approximation initiale et de f
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
elseif ~isnumeric(x1) || ~isscalar(x1)
	error('L''approximation initiale x1 n''est pas un scalaire')
end

try
	fct(x0);
catch ME
	rethrow(ME)
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0))
	error('Le fonction f ne retourne pas un scalaire')
end

%% Vérification du fichier output
if nargin == 6 && ~isa(file_name,'char')
	error('Le nom du fichier des résultats doit être de type string')
end

%% Initialisation des matrices app et err
app			=	nan(1,nb_it_max);
app(1)		=	x0;
app(2)		=	x1;
err_rel		=	inf(1,nb_it_max);
err_rel(1)  =	(x1 - x0)/(x1+eps);
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
err_abs		=	inf(1,nb_it);

if arret
	err_abs		=	abs(approx(end) - approx);
else
	warning('La méthode de la sécante n''a pas convergée')
end

% Écriture des résultats si fichier passé en argument
if nargin == 6
	output_results(file_name , fct , is_fct_file , nb_it_max , ...
					tol_rel , x0 , x1 , approx , err_abs , arret)
end

end

function [] = output_results(file_name , fct , is_fct_file , it_max , ...
							tol_rel , x0 , x1 , x , err , status)
						 
	fid		=	fopen(file_name,'w');
	fprintf(fid,'Algorithme de la sécante\n\n');
	fprintf(fid,'Fonction dont on cherche les racines:\n');
	if is_fct_file
		fprintf(fid,'%s\n\n',fileread([func2str(fct),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(fct));
	end
	fprintf(fid,'Arguments d''entrée:\n');
	fprintf(fid,'    - Nombre maximum d''itérations: %d\n',it_max);
	fprintf(fid,'    - Tolérance relative: %6.5e\n',tol_rel);
	fprintf(fid,'    - Approximation initiale x0: %16.15e\n',x0);
	fprintf(fid,'    - Approximation initiale x1: %16.15e\n\n',x1);
	
	if status
		fprintf(fid,'\nStatut: L''algorithme de la sécante a convergé en %d itérations\n\n',length(x)-2);
	else
		fprintf(fid,'\nStatut: L''algorithme de la sécante n''a pas convergé\n\n');
	end
	fprintf(fid,'#It            x              Erreur absolue\n');
	fprintf(fid,'%3d   %16.15e   %6.5e\n',[reshape(0:length(x)-1,1,[]);reshape(x,1,[]);reshape(err,1,[])]);
	fclose(fid);

end
