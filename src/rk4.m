function [temps , y] = rk4(f , tspan , Y0 , nb_pas)
% RK4	Méthode de Runge-Kutta d'ordre 4 avec pas constant pour la
%		résolution d'EDOs
%
% Syntaxe: [temps , y] = rk4(f , tspan , Y0 , nb_pas)
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
%	[temps , y] = rk4(@(t,y) y*cos(t) , [0,2] , 1 , 1000 );
%	[temps , y] = rk4(@(t,z) [z(2);-10*z(1)] , [0,1] , [1;0] , 1000);
%	[temps , y] = rk4('my_edo' , [0,1] , [1;0] , 1000);



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
if ~isnumeric(Y0) || ~isvector(Y0)
	error('Les conditions initiales x0 ne sont pas arrangées en vecteur')
end

try 
	fct(t0,Y0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction f n''est pas dans le répertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le nombre de composantes de x0 et f ne concorde pas')
	else 
		rethrow(ME)
	end
end

if ~isnumeric(fct(t0,Y0)) || ~isvector(fct(t0,Y0))
	error('La f ne retourne pas un vecteur')
elseif length(Y0) ~= length(fct(t0,Y0))
	error('Le nombre de composantes de x0 et f ne concorde pas')
end

nb_comp		=	length(Y0);

%% Initialisation du temps et de la matrice y
temps		=	reshape(linspace(t0,tf,nb_pas+1),nb_pas+1,1);
y			=	nan(nb_pas + 1 , length(Y0));
y(1,:)		=	reshape(Y0,1,nb_comp);


%% Méthode de Runge-Kutta d'ordre 4 avec pas constant
for t=1:nb_pas
	k1		=	h*reshape(fct(temps(t)		 , y(t,:))		  ,1,nb_comp);
	k2		=	h*reshape(fct(temps(t) + h/2 , y(t,:) + k1/2 ),1,nb_comp);
	k3		=	h*reshape(fct(temps(t) + h/2 , y(t,:) + k2/2 ),1,nb_comp);
	k4		=	h*reshape(fct(temps(t) + h   , y(t,:) + k3   ),1,nb_comp);
	y(t+1,:)=	y(t,:) + 1/6 * ( k1 + 2*k2 + 2*k3 + k4 );
end

end