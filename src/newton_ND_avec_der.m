function [approx , err_abs] = newton_ND_avec_der(F , mat_jac , x0 , nb_it_max , tol_rel , file_name)
% NEWTON_ND_AVEC_DER	Méthode de Newton pour la résolution de F(x) = 0, pour F: R^n -> R^n
%
% Syntaxe: [approx , err_abs] = newton_ND_avec_der(f , Jac , x0 , nb_it_max , tol_rel)
%
% Argument d'entrée
%	F			-	String ou fonction handle spécifiant la fonction
%					non-linéaire (F: R^n -> R^n)
%	mat_jac		-	String ou function handle spécifiant la matrice
%					jacobienne de F (Jac: R^n -> R^{n x n})
%	x0			-	Approximation initiale (x0 in R^n)
%	nb_it_max	-	Nombre maximum d'itérations 
%	tol			-	Tolérance sur l'approximation de l'erreur relative
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
%	[ approx , err_abs ] = newton_ND_avec_der( 'my_sys_nl' , 'my_sys_nl_jac' , [1,1] , 20 , 1e-9 , 'resul_newtonND.txt')
%	[ approx , err_abs ] = newton_ND_avec_der( @(x) [x(1)^2 + x(2)^2 - 1 ; -x(1)^2 + x(2)] , @(x) [2*x(1),2*x(2);-2*x(1),1] , [1,1] , 20 , 1e-9 , 'resul_newtonND.txt')



%%  Vérification de la fonction contenant les dérivées
if isa(F,'char')
	fct			=	str2func(F);
	is_fct_file =	true;
elseif isa(F,'function_handle')
	fct			=	F;
	is_fct_file =	false;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

if isa(mat_jac,'char')
	jac				=	str2func(mat_jac);
	is_jac_file		=	true;
elseif isa(mat_jac,'function_handle')
	jac				=	mat_jac;
	is_jac_file		=	false;
else
	error('L''argument df n''est pas un string ni un function_handle')
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

try 
	jac(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction mat_jac n''est pas dans le répertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction mat_jac ne retourne pas une matrice carré')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction mat_jac doit prendre seulement 1 argument en entrée')
	else
		rethrow(ME)
	end
end

taille	=	length(x0);

if ~isnumeric(fct(x0)) || ~isvector(fct(x0)) || (length(fct(x0))~=taille)
	error(['Le fonction F ne retourne pas un vecteur de même taille que x0. ',...
		'x0 est de taille %d alors que F(x0) est de taille %d.'],taille,length(fct(x0)))
elseif ~isnumeric(jac(x0)) || ~ismatrix(jac(x0)) || ~isequal(size(jac(x0)),[taille,taille]) 
	error(['Le fonction mat_jac ne retourne pas une matrice de bonne taille; ',...
		   'x0 est de taille %d, alors que mat_jac(x0) est de dimension %d x %d.'],...
										taille,size(jac(x0),1),size(jac(x0),2))
elseif ~check_derivative(fct,jac,x0)
	warning('Il semble y avoir une erreur avec la matrice jacobienne')
end

%% Initialisation des matrices app et err
app			=	nan(nb_it_max,taille);
app(1,:)	=	x0;
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% Méthode de Newton
for t=1:nb_it_max-1
	
	delta_x		=	jac(app(t,:))\-reshape(fct(app(t,:)),taille,1);
	app(t+1,:)	=	app(t,:) + delta_x';
	
	if min(abs(eig(jac(app(t,:))))) == 0  
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

if nargin == 6
	output_results(file_name , fct , is_fct_file , jac, is_jac_file, ...
				nb_it_max , tol_rel , x0 , approx , err_abs , arret)
end

end


function [test] = check_derivative(f,jac,x0)

	taille	=	length(x0);
	if min(x0) == 0
		h_init	=	1e-6;
	else
		h_init	=	1e-3 * min(x0);
	end	
	h		=	h_init./(2.^(0:4));
	erreur	=	nan(size(h));

	
	for t=1:length(h)
		app			=	zeros(taille,taille);
		for d=1:taille
			delta_h		=	zeros(size(x0));
			delta_h(d)	=	h(t);
			app(:,d)	=	(f(x0+delta_h) - f(x0-delta_h))/(2*h(t));
		end
		erreur(t)	=	norm(jac(x0) - app);
	end
	
	ordre	=	log(erreur(2:end) ./ erreur(1:end-1)) ./ ...
				log(h(2:end)      ./ h(1:end-1));
			
	test	=	(mean( abs(2 - ordre) ) <= 0.25) || (mean(erreur/(norm(jac(x0))+eps)) <=1e-9);	
	
end

function [] = output_results(file_name , fct , is_fct_file , jac, ...
			is_jac_file, it_max , tol_rel , x0 , x , err , status)
						 
	fid		=	fopen(file_name,'w');
	fprintf(fid,'Algorithme de Newton avec la matrice jacobienne\n\n');
	fprintf(fid,'Fonction dont on cherche les racines:\n');
	if is_fct_file
		fprintf(fid,'%s\n\n',fileread([func2str(fct),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(fct));
	end
	fprintf(fid,'Matrice jacobienne de la fonction dont on cherche les racines:\n');
	if is_jac_file
		fprintf(fid,'%s\n\n',fileread([func2str(jac),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(jac));
	end
	fprintf(fid,'Arguments d''entrée:\n');
	fprintf(fid,'    - Nombre maximum d''itérations: %d\n',it_max);
	fprintf(fid,'    - Tolérance relative: %6.5e\n',tol_rel);
	fprintf(fid,'    - Approximation initiale x0: [');
	taille	=	size(x,2);
	if taille <= 3
		fprintf(fid,'%6.5e ',x0(1:taille));
		fprintf(fid,']\n');
	else
		fprintf(fid,'%6.5e ',x0(1:2));
		fprintf(fid,'... %6.5e]',x0(end));
	end
	
	if status
		fprintf(fid,'\nStatut: L''algorithme de Newton a convergé en %d itérations\n\n',size(x,1)-1);
	else
		fprintf(fid,'\nStatut: L''algorithme de Newton n''a pas convergé\n\n');
	end
	if taille == 1
		fprintf(fid,'#It       x_1           Erreur absolue\n');
		fprintf(fid,'%3d   %16.15e   %6.5e\n',[reshape(0:length(x)-1,1,[]);reshape(x(:),1,[]);reshape(err,1,[])]);
	elseif taille == 2
		fprintf(fid,'#It       x_1           x_2       Erreur absolue\n');
		fprintf(fid,'%3d   %6.5e   %6.5e   %6.5e\n',[reshape(0:length(x)-1,1,[]);x(:,[1,2])';reshape(err,1,[])]);
	elseif taille == 3
		fprintf(fid,'#It       x_1           x_2           x_3       Erreur absolue\n');
		fprintf(fid,'%3d   %6.5e   %6.5e   %6.5e\n',[reshape(0:length(x)-1,1,[]);x(:,[1,2,3])';reshape(err,1,[])]);
	else
		fprintf(fid,'#It       x_1           x_2      ...       x_n       Erreur absolue\n');
		fprintf(fid,'%3d   %6.5e   %6.5e  ...  %6.5e\n',[reshape(0:length(x)-1,1,[]);x(:,[1,2,taille])';reshape(err,1,[])]);
	end
	
	fclose(fid);

end
