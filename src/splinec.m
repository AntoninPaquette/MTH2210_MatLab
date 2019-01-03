function [ Sx ] = splinec( xi , yi , x , type_S , val_S)
% SPLINEC	Spline cubique passant par les points xi et yi avec différents 
%			type de conditions frontières
%
% Syntaxe: [ Sx ] = splinec( xi , yi , x , type_f , val_f)
%
% Argument d'entrée
%	xi		-	Vecteur contenant les abscisses des points d'interpolation
%	yi		-	Vecteur contenant les ordonnées des points d'interpolation
%	x		-	Vecteur contenant les points où la spline sera évalué
%	type_S	-	Vecteur de dimension 2 contenant le type des conditions 
%				frontières imposées en x0 et xn. Les choix possibles sont:
%					[1,1] -> Spline naturelle 
%					[2,2] -> Spline avec courbure prescrite
%					[3,3] -> Spline avec courbure constante
%					[4,4] -> Spline avec pente prescrite
%					[i,j] -> Spline avec condition i imposée en x0 et 
%							 condition j imposée en xn
%							 	
%	val_S	-	Vecteur de dimension 2 contenant les deux conditions 
%				limites imposées en x0 et xn. Les choix possibles sont:
%					- Si type_S(1) = 1 ou 3, alors val_S(1) = nan
%					- Si type_S(1) = 2 ou 4, alors val_S(1) = a, où a 
%					  représente resp. la courbure ou la pente en x0
%					- Si type_S(2) = 1 ou 3, alors val_S(1) = nan
%					- Si type_S(2) = 2 ou 4, alors val_S(1) = b, où b 
%					  représente resp. la courbure ou la pente en xn
%
% Arguments de sortie
%	Sx		-	Vecteur contenant les valeurs de la spline aux points x
%
% Exemples d'appel
%	[ Sx ] = splinec([1,2,4,5], [1,9,2,11], linspace(1,5), [1,1] , [nan,nan])
%	[ Sx ] = splinec([1,2,4,5], [1,9,2,11], linspace(1,5), [2,2] , [5,-6])
%	[ Sx ] = splinec([1,2,4,5], [1,9,2,11], linspace(1,5), [3,3] , [nan,nan])
%	[ Sx ] = splinec([1,2,4,5], [1,9,2,11], linspace(1,5), [4,4] , [-30,-10])
%	[ Sx ] = splinec([1,2,4,5], [1,9,2,11], linspace(1,5), [3,4] , [nan,-10])


%% Vérification des arguments d'entrées
if ~isnumeric(xi) || ~isvector(xi)
	error('Les abscisses xi ne sont pas arrangées en vecteur')
elseif ~isequal(sort(xi),xi)
	error('Les abscisses ne sont pas arrangées en ordre croissant')
elseif ~isnumeric(yi) || ~isvector(yi)
	error('Les ordonnées yi ne sont pas arrangées en vecteur')
elseif length(xi) ~= length(yi)
	error('Les vecteurs xi et yi doivent avoir la meme taille');
elseif ~isnumeric(x) || ~isvector(x)
	error('Les points où l''on évalue la spline cubique doivent être arrangées en vecteur');
elseif isequal(xi,x)
	warning('Le polynôme d''interpolation est évalué exactement au points d''interpolation')
elseif ~isnumeric(type_S) || ~isvector(type_S) || length(type_S)~=2 || ~any(type_S(1)==[1,2,3,4]) || ~any(type_S(2)==[1,2,3,4])
	error('Les types de conditions frontières ne sont pas valide')
elseif ~isnumeric(val_S) || ~isvector(val_S) || length(val_S)~=2
	error('Les valeurs des conditions frontières ne sont pas dans un vecteur de dimension 2')
end


%% Calcul des coefficients S''

% Assemblage de la matrice
nb_f	=	length(xi);
h		=	reshape(diff(xi),1,[]);
denom	=	h(1:end-1)+h(2:end);

M			=	[[h(1:end-1)./denom,0,0]',[0;2*ones(nb_f-2,1);0],[0,0,h(2:end)./denom]'];
matrice		=	spdiags(M,[-1,0,1],nb_f,nb_f);

% Calcul des 2ème différences divisées
mat_diff_div1	=	diff_div(xi,yi,1);
mat_diff_div2	=	diff_div(xi,mat_diff_div1,2);

% Terme de droite
B	=	zeros(nb_f,1);
B(2:end-1)	=	6*mat_diff_div2;

% Imposition des conditions frontières
switch type_S(1)
	case 1
		matrice(1,1)		=	1;
		B(1)				=	0;
	case 2
		matrice(1,1)		=	1;
		B(1)				=	val_S(1);
	case 3
		matrice(1,[1,2])		=	[1,-1];
		B(1)					=	0;
	case 4
		matrice(1,[1,2])		=	[2,1];
		B(1)					=	6/h(1)   * ((yi(2) - yi(1))/h(1) - val_S(1));
end

switch type_S(2)
	case 1
		matrice(nb_f,nb_f)	=	1;
		B(end)				=	0;
	case 2
		matrice(nb_f,nb_f)	=	1;
		B(end)				=	val_S(2);
	case 3
		matrice(end,[end-1,end])=	[-1,1];
		B(end)					=	0;
	case 4
		matrice(end,[end-1,end])=	[1,2];
		B(end)					=	6/h(end) * (val_S(2) - (yi(end) - yi(end-1))/h(end));
end

% Résolution du système linéaire
Spp		=	matrice\B;


%% Calcul de la spline aux points x 

Sx	=	nan(size(x));

for t=1:nb_f-1
	x_inter		=	(x>= xi(t)) & (x <= xi(t+1));
	Sx(x_inter) =	poly_spline(xi(t:t+1),yi(t:t+1),Spp(t:t+1),h(t),x(x_inter));
end

end


function [ f ] = diff_div(x,y,n)
	diff_x	=	reshape(x(n+1:end) - x(1:end-n),1,[]);
	diff_y	=	reshape(diff(y),1,[]);
	f		=	diff_y./diff_x;
end


function [ Px ] = poly_spline(xi,yi,Spp,h,x)
	Px	=	-Spp(1)/(6*h)*(x-xi(2)).^3 + Spp(2)/(6*h)*(x-xi(1)).^3 - ...
			 (yi(1)/h - Spp(1)*h/6)*(x-xi(2)) + (yi(2)/h - Spp(2)*h/6)*(x-xi(1));
end