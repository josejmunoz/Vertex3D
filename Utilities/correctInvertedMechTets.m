function [Y, Cell] = correctInvertedMechTets(exceptionMessage, dy, Y, Cell, Set)
%CORRECTINVERTEDMECHTETS Corrects inverted 'mechanical' tetrahedra (don't
%confused with 'geometrical' tetrahedra T)
%   Detailed explanation goes here

disp('Correcting inverted mechanical tetrahedra...')

indicesToUpdate = str2num(exceptionMessage.message(31:end));
        
% Assemble numbers
dim = 3;
idofg = zeros(length(indicesToUpdate)*dim, 1);

for I=1:length(indicesToUpdate) % loop on number of vertices of triangle
    idofg((I-1)*dim+1:I*dim) = (indicesToUpdate(I)-1)*dim+1:indicesToUpdate(I)*dim;% global dof
end

dy_reajusted = zeros(size(dy));
dy_reajusted(idofg) = - dy(idofg);

dy_reshaped = reshape(dy_reajusted, 3, Set.NumTotalV)';

[Y, Cell] = updateVertices(Y, Cell, dy_reshaped, Set);
end

