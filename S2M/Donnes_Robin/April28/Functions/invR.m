function [R10] = invR(R01, TwoDimension)
%% Auteur :                 BEGON Mickael
%% Organisation :           Laboratoire de M�canique des Solides, Poitiers
%% Cr�ation :               

% simple fonction pour passer de la matrice de passage R01 � son inverse R10

% % R = R01(1:3,1:3);
% % T = R01(1:3,4);
% % 
% % R10 = [R01'; 0 0 0];
% % R10 = [R01 [-R01*T;1]];

    if isnumeric(R01)
        [n,~,q] = size(R01);
        
        if nargin >=2 && TwoDimension
            switch n
                case 3
                    R = multitransp(R01(1:2,1:2,:));
                    T = R01(1:2,3,:);
                    R10 = [R, multiprod(-R,T); [zeros(1,2,q) ones(1,1,q)]];
                case 2
                    R10 = multitransp(R01);
            end
            
        else
            switch n
                case 4
                    R = multitransp(R01(1:3,1:3,:));
                    T = R01(1:3,4,:);
                    R10 = [R, multiprod(-R,T); [zeros(1,3,q) ones(1,1,q)]];
                case 3
                    R10 = multitransp(R01);
            end
        end
    else % calcul symbolique
        R = [R01(1,1) R01(2,1) R01(3,1); R01(1,2) R01(2,2) R01(3,2); R01(1,3) R01(2,3) R01(3,3)];
        T = R01(1:3,4);
        R10 = [R -R*T; 0 0 0 1];
    end
                
end

