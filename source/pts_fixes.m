function [approx , err_abs] = pts_fixes(g , x0 , nb_it_max , tol_rel , file_name)
% PTS_FIXES	Methode des points-fixes pour la resolution de g(x) = x
%			pour g: R -> R
%
% Syntaxe: [approx , err_abs] = pts_fixes(g , x0 , nb_it_max , tol_rel , file_name)
%
% Arguments d'entree
%	g			-	String ou fonction handle specifiant la fonction
%					non-lineaire
%	x0			-	Approximation initiale 
%	nb_it_max	-	Nombre maximum d'iterations 
%	tol_rel		-	Tolerance sur l'approximation de l'erreur relative
%	file_name	-	(Optionnel) Nom du fichier (avec l'extension .txt) dans
%					lequel sera	ecrit les resultats de l'algorithme
%
% Arguments de sortie
%	approx		-	Vecteur rangee de taille nb_iter contenant les 
%					iterations
%	err_abs		-	Vecteur rangee de dimension nb_iter contenant les
%					erreurs absolues
%
% Exemples d'appel
%	[approx , err_abs] = pts_fixes('my_fct_ptfixe' , 4 , 50 , 1e-12 , 'resul_pts_fixes.txt')
%	[approx , err_abs] = pts_fixes(@(x) -x.^2/10+x+1 , 4 , 50 , 1e-12 , 'resul_pts_fixes.txt')



%%  Verification de la fonction g
if isa(g,'char')
	fct			=	str2func(g);
	is_fct_file =	true;
elseif isa(g,'function_handle')
	fct			=	g;
	is_fct_file =	false;
else
	error('L''argument g n''est pas un string ni un function_handle')
end

%% Verification du nb de composantes de l'approximation initiale et de g
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
end

try 
	fct(x0);
catch ME
	rethrow(ME)
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0)) 
	error('Le fonction g ne retourne pas un scalaire')
end

%% Verification du fichier output
if nargin == 5 && ~isa(file_name,'char')
	error('Le nom du fichier des resultats doit etre de type string')
end

%% Initialisation des matrices app et err
app			=	nan(1,nb_it_max);
app(1)		=	x0;
err_rel		=	inf(1,nb_it_max);
arret		=	false;


%% Methode des points fixes
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
err_abs		=	inf(1,nb_it);

if arret
	err_abs		=	abs(approx(end) - approx);
else
	warning('La methode des points-fixes n''a pas convergee')
end

% ecriture des resultats si fichier passe en argument
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

	fprintf(fid,'Arguments d''entree:\n');
	fprintf(fid,'    - Nombre maximum d''iterations: %d\n',it_max);
	fprintf(fid,'    - Tolerance relative: %6.5e\n',tol_rel);
	fprintf(fid,'    - Approximation initiale x0: %16.15e\n',x0);
	
	if status
		fprintf(fid,'\nStatut: L''algorithme des points-fixes a converge en %d iterations\n\n',length(x)-1);
	else
		fprintf(fid,'\nStatut: L''algorithme des points-fixes n''a pas converge\n\n');
	end
	fprintf(fid,'#It            x              Erreur absolue\n');
	fprintf(fid,'%3d   %16.15e   %6.5e\n',[reshape(0:length(x)-1,1,[]);reshape(x,1,[]);reshape(err,1,[])]);
	fclose(fid);

end