function [] = affichage_tableau(varargin)

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
            error("Les vecteurs doivent êtres de même taille")
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
    fprintf(titre);
    fprintf(join(repmat("-",1,taille_tot),"")+"\n")
    
    % Affichage des données
    for rows=1:nb_rows
        fprintf("| ")
        for cols=1:nb_cols
            if is_dec(cols)
                fprintf("%*d | ",taille(cols),varargin{2*cols-1}(rows))
            else
                fprintf("%22.15e | ",varargin{2*cols-1}(rows));
            end    
        end
        fprintf("\n");
    end
    
end

