function [psi_out, theta_out, phi_out] = angleRotation(matrice, sequence)
%Cette fonction retourne les angles psi, th�ta et phi d'une matrice de rotation en fonction 
%de la s�quence d'angle souhait�e. 
%�crite par pariterre le 15 f�vrier 2011
%Il est � not� que le tout a �t� trouv� � la main. J'ai contr�l� au maximum pour les erreurs, mais 
%il est toujours possible qu'il en reste
% [psi, theta, phi] : DOUBLE, contenant les angles. Si 1 seul argument de sorti : tout est plac� dans psi
% matrice : MATRICE DOUBLE, 3*3 
% sequence : STRING, de la s�quence d'angle


    %Controle du nombre d'argument
    if nargin ~= 2 
        error('Arguments manquants dans la fonction angleRotation')
    end %if
    
    if size(matrice,4) > 1 && (size(matrice,4) ~= length(sequence))
        error('Le nombre de séquence doit correspondre au nombre d''éléments à calculer')
    elseif size(matrice,4) == 1
        sequence = {sequence}; % Transformer en cellule
    end
        
    psi_out = nan(1,size(matrice,4), size(matrice,3));
    theta_out = nan(1,size(matrice,4), size(matrice,3));
    phi_out = nan(1,size(matrice,4), size(matrice,3));
    for iM = 1:size(matrice, 4)
        psi = []; %#ok<NASGU>
        theta = [];
        phi = [];
        switch sequence{iM}
            case 'x'
                psi = asin(matrice(3,2,:,iM));
            case 'y'
                psi = asin(matrice(1,3,:,iM));
            case 'z'
                psi = asin(matrice(2,1,:,iM));
            case 'xy'
                psi = asin(matrice(3,2,:,iM));
                theta = asin(matrice(1,3,:,iM));
            case 'xz'
                psi = -asin(matrice(2,3,:,iM));
                theta = -asin(matrice(1,2,:,iM));
            case 'yx'
                psi = -asin(matrice(3,1,:,iM));
                theta = -asin(matrice(2,3,:,iM));
            case 'yz'
                psi = asin(matrice(1,3,:,iM));
                theta = asin(matrice(2,1,:,iM));
            case 'zx'
                psi = asin(matrice(2,1,:,iM));
                theta = asin(matrice(3,2,:,iM));
            case 'zy'
                psi = -asin(matrice(1,2,:,iM));
                theta = -asin(matrice(3,1,:,iM));
            case 'xyz'
                psi   = atan2(-matrice(2,3,:,iM), matrice(3,3,:,iM)); 	
                theta = asin(matrice(1,3,:,iM)); 				
                phi   = atan2(-matrice(1,2,:,iM), matrice(1,1,:,iM));
            case 'xzy'   
                psi   = atan2(matrice(3,2,:,iM),matrice(2,2,:,iM));       
                phi = atan2(matrice(1,3,:,iM), matrice(1,1,:,iM));   
                theta   = asin(-matrice(1,2,:,iM));
            case 'yxz'
                theta   = asin(-matrice(2,3,:,iM));
                psi = atan2(matrice(1,3,:,iM), matrice(3,3,:,iM));
                phi   = atan2(matrice(2,1,:,iM), matrice(2,2,:,iM));
            case 'yzx'
                phi   = atan2(-matrice(2,3,:,iM), matrice(2,2,:,iM));
                psi = atan2(-matrice(3,1,:,iM), matrice(1,1,:,iM));
                theta   = asin(matrice(2,1,:,iM));        
            case 'zxy'
                theta   = asin(matrice(3,2,:,iM));
                phi = atan2(-matrice(3,1,:,iM), matrice(3,3,:,iM));
                psi   = atan2(-matrice(1,2,:,iM), matrice(2,2,:,iM));
            case 'zyx'
                phi   = atan2(matrice(3,2,:,iM), matrice(3,3,:,iM));
                theta = asin(-matrice(3,1,:,iM));
                psi   = atan2(matrice(2,1,:,iM), matrice(1,1,:,iM));
            case 'zyz'
                psi   = atan2(matrice(2,3,:,iM), matrice(1,3,:,iM));
                theta = acos(matrice(3,3,:,iM));
                phi   = atan2(matrice(3,2,:,iM), -matrice(3,1,:,iM));
            case 'zxz'
                psi	  = atan2(matrice(1,3,:,iM), -matrice(2,3,:,iM));
                theta = acos(matrice(3,3,:,iM));
                phi   = atan2(matrice(3,1,:,iM), matrice(3,2,:,iM));
            case 'zyzz'
                psi   = atan2(matrice(2,3,:,iM), matrice(1,3,:,iM));
                theta = acos(matrice(3,3,:,iM));
                phi   = atan2(matrice(3,2,:,iM), -matrice(3,1,:,iM)) + psi;


            otherwise
                error('Séquence d''angle incorrecte dans angleRotation')

        end
        if ~isempty(psi)
            psi_out(:,iM,:) = real(psi);
        end
        if ~isempty(theta)
            theta_out(:,iM,:) = real(theta);
        end
        if ~isempty(phi)
            phi_out(:,iM,:) = real(phi);
        end
    end

    if nargout == 1			%Si un seul argument a �t� envoy�, mettre tous les angles dans une seule variable
        psi_out = [psi_out; theta_out; phi_out];
        if length(sequence) == 1 % Retitrer les nan que s'il y a un segment
            if strcmp(sequence{1}, 'zyzz')
                psi_out = psi_out(1:3,:,:);
            else
                psi_out = psi_out(1:length(sequence{1}),:,:);
            end
        end
    end
end %function
