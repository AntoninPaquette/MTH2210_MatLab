function [ Lx ] = lagrange(xi , yi , x)
% LAGRANGE	Polynôme de Lagrange passant par les points xi et yi
%
% Syntaxe: [ Ly ] = lagrange(xi , yi , x)
%
% Argument d'entrée
%	xi		-	Vecteur contenant les abscisses des points d'interpolation
%	yi		-	Vecteur contenant les ordonnées des points d'interpolation
%	x		-	Vecteur contenant les points où le polynôme de Lagrange 
%				sera évalué
%
% Arguments de sortie
%	Lx		-	Vecteur contenant les valeurs du polynômes aux points x
%
% Exemples d'appel
%	[ Lx ] = lagrange([-1,0,1] , [1,0,1] , linspace(-1,1,200))
%
% Le calcul du polynome interpolant est base sur la formule barycentrique ;
% une variation stable de la formule d'interpolation de Lagrange vue au
% cours. Cette implementation est basee sur celle de Greg von Winckel,
% elle-meme basee sur l'article de Berrut et Trefethen (SIAM Review, 2004).
%
% Normalement, on devrait utiliser une version particularisee de cette
% fonction si le vecteur xi donne les noeuds d'interpolation de Chebychev.
%
% "Barycentric interpolation is a variant of Lagrange polynomial
% interpolation that is fast and stable. It deserves to be known as the
% standard method of polynomial interpolation." (Berrut and Trefethen, 2004).
%
% Créé par D. Orban, 2009. 
% Modifié par A. Paquette-Rufiange, 2018.
%

%% Vérification des arguments d'entrées
if ~isnumeric(xi) || ~isvector(xi)
	error('Les abscisses xi ne sont pas arrangées en vecteur')
elseif ~isnumeric(yi) || ~isvector(yi)
	error('Les ordonnées yi ne sont pas arrangées en vecteur')
elseif length(xi) ~= length(yi)
	error('Les vecteurs xi et yi doivent avoir la meme taille');
elseif ~isnumeric(x) || ~isvector(x)
	error('Les points où l''on évalue la spline cubique doivent être arrangées en vecteur');
elseif isequal(xi,x)
	warning('Le polynôme d''interpolation est évalué exactement au points d''interpolation')
end


%% Calcul de l'interpolant aux points x

M = length(xi);  xxi = reshape(xi,M,1);  % Donnees en vecteur colonne
N = length(x);   xx = reshape(x,N,1);    % Grille fine en vecteur colonne

Xi = repmat(xxi,1,M);
W = repmat(1./prod(Xi-Xi.'+eye(M),1),N,1);   % Poids barycentriques
xdist=repmat(xx,1,M)-repmat(xxi.',N,1);   % Distances aux points d'interpolation

% On place des NaN aux noeuds d'interpolation pour eviter la division par zero
[fixi,fixj]=find(xdist==0);
xdist(fixi,fixj)=NaN;
H=W./xdist;

Lx=(H * reshape(yi,M,1))./sum(H,2);  % Valeurs du polynome sur la grille fine
Lx(fixi)=yi(fixj);        % On restaure les valeurs exacte a la place des NaN
if (size(Lx) ~= size(x))
  Lx = Lx';
end
