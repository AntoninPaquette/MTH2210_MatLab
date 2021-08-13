function [str] = affichage_tableau(varargin)
% AFFICHAGE_TABLEAU	Fonction permettant de créer un string mettant sous
%					format un tableau de valeur avec leur titre
%
% Syntaxe: [str] = affichage_tableau(vecteur_1,titre_1,vecteur_2,titre_2)
%
% Arguments d'entrée  
%	vecteur_i	-	Vecteur de donnees a afficher à la i-eme colonne 
% 					(premier argument d'une paire) 
%	titre_i		-	String spécifiant le titre de la i-eme colonne
%					(deuxieme argument d'une paire)
%
% Arguments de sortie
%	str			-	String contenant le tableau pouvant être affiche avec
%					la fonction fprintf 
%
% Exemples d'appel
%	[ str ] = affichage_tableau(0:2:10,"n",10.^(0:2:10),"10^n")
%	[ str ] = affichage_tableau(linspace(0,5,15),"x",cos(linspace(0,5,15)),"cos(x)")


% Vérification du nombre d'inputs
if rem(nargin,2) ~= 0
	error("Les entrées doivent être en paire vecteur/titre.")
end
nb_cols = nargin/2;

% Vérification des inputs
for t=1:nb_cols
	if ~isvector(varargin{2*t-1})
		error("Les premiers éléments des paires doivent être des vecteurs.")
	end
	if ~isa(varargin{2*t},"string")
		error("Les deuxièmes éléments des paires doivent être des strings.")
	end
end

% Vérification de la taille des vecteurs
nb_rows = length(varargin{1});
for t=1:nb_cols
	if length(varargin{2*t-1}) ~= nb_rows
		error("Les vecteurs doivent être de même taille")
	end
end

% Formattage et titre des colonnes
is_dec = false(nb_cols,1);
taille = 22*ones(nb_cols,1);
titre = "| ";
for t=1:nb_cols
	is_dec(t) = all(mod(varargin{2*t-1},1)==0);
	if is_dec(t)
		taille_string = strlength(varargin{2*t});
		if any(varargin{2*t-1}<0)
			taille_dec = max(floor(log10(abs(varargin{2*t-1})))) + 2;
		else
			taille_dec = max(floor(log10(abs(varargin{2*t-1})))) + 1;
		end
		taille(t) = max(taille_string,taille_dec);
		titre = titre + sprintf("%-*s | ",taille(t),varargin{2*t});
	else
		titre = titre + sprintf("%-22s | ",varargin{2*t});
	end
end
titre = titre + "\n";

% Vérification de la taille totale du tableau
taille_tot = sum(taille) + 3*nb_cols + 1;
if taille_tot > 109
	warning("Le tableau affiché est peut-être trop large.")
end

% Affichage du titre
str		=	sprintf(titre);
str		=	str + sprintf(join(repmat("-",1,taille_tot),"")+"\n");

% Affichage des données
for rows=1:nb_rows
	str		=	str + sprintf("| ");
	for cols=1:nb_cols
		if is_dec(cols)
			str		=	str + sprintf("%*d | ",taille(cols),varargin{2*cols-1}(rows));
		else
			str		=	str + sprintf("%22.15e | ",varargin{2*cols-1}(rows));
		end
	end
	str		=	str + newline;
end

end

