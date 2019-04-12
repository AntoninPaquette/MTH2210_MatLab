function [temps , y] = crank_nic(f , tspan , x0 , nb_pas)
% CRANK_NIC	Méthode de Crank-Nicolson avec pas constant pour la
%			résolution d'EDOs
%
% Syntaxe: [temps , y] = crank_nic(f , tspan , Y0 , nb_pas)
%
% Arguments d'entrée
%	f		-	String ou function handle définissant le système de N EDOs
%	tspan	-	Vecteur contenant le temps initial et final [t0,tf]
%	x0		-	Vecteur contenant les N conditions initiales
%	nb_pas	-	Nombre de pas de temps
%
% Arguments de sortie
%	temps	-	Vecteur colonne contenant les valeurs de temps t_i
%	y		-	Matrice de dimension (nb_pas+1) x N dont les colonnes 
%				sont les approximations de y_i(t)
%
% Exemples d'appel
%	[temps , y] = crank_nic('my_edo' , [0,1] , [1;0] , 1000);
%	[temps , y] = crank_nic(@(t,y) y*cos(t) , [0,2] , 1 , 1000 );
%	[temps , y] = crank_nic(@(t,z) [z(2);-10*z(1)] , [0,1] , [1;0] , 1000);


%%  Vérification de la fonction contenant les dérivées
if isa(f,'char')
	fct		=	str2func(f);
elseif isa(f,'function_handle')
	fct		=	f;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

%% Vérification du temps et nb pas de temps
if ~isnumeric(tspan) || length(tspan)~=2
	error('Le vecteur tspan doit contenir 2 composantes, [t0 , tf]')
elseif ~isnumeric(nb_pas) || floor(nb_pas)~=nb_pas ...
						  || length(nb_pas)~=1 || nb_pas<0
	error('Le nombre de pas nb_pas doit être entier et positif')
end

t0			=	tspan(1);
tf			=	tspan(2);
nb_pas		=	double(nb_pas);
h			=	(tf-t0)/nb_pas;

%% Vérification du nb de composantes des conditions initiales et de f
if ~isnumeric(x0) || ~isvector(x0)
	error('Les conditions initiales x0 ne sont pas arrangées en vecteur')
end

try 
	fct(t0,x0);
catch ME
	rethrow(ME)
end

if ~isnumeric(fct(t0,x0)) || ~isvector(fct(t0,x0))
	error('La f ne retourne pas un vecteur')
elseif length(x0) ~= length(fct(t0,x0))
	error('Le nombre de composantes de x0 et f ne concorde pas')
end

nb_comp		=	length(x0);

%% Initialisation du temps et de la matrice y
temps		=	reshape(linspace(t0,tf,nb_pas+1),nb_pas+1,1);
y			=	nan(nb_pas + 1 , length(x0));
y(1,:)		=	reshape(x0,1,nb_comp);

%% Méthode d'Euler implicite avec pas constant
for t=1:nb_pas
	fct_nl	=	@(x) reshape(x,1,nb_comp) - y(t,:) - h/2*reshape(fct(temps(t),y(t,:)) + fct(temps(t+1),x),1,nb_comp); 
	[approx , err_abs] = newton_ND_sans_der(fct_nl , y(t,:) , 10 , 1e-9);
	if isinf(err_abs(end))
		warning('Problème au temps %1.6e',temps(t+1))
	end
	y(t+1,:)	=	approx(end,:);
end


end

