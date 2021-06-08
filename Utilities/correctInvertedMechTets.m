function [] = correctInvertedMechTets(Cell, Y)
%CORRECTINVERTEDMECHTETS Corrects inverted 'mechanical' tetrahedra (don't
%confused with 'geometrical' tetrahedra T)
%   Detailed explanation goes here


%% K and g calculation per Cell
for numCell = 1:ncell
    cellNuclei = Cell.Centre(numCell, :);
    % Loop over Cell-face-triangles
    Tris=Cell.Tris{numCell};
    for ntriangle=1:size(Tris,1)
        currentTet_ids=Tris(ntriangle,:);

        Y1=Y.DataRow(currentTet_ids(1),:);
        Y2=Y.DataRow(currentTet_ids(2),:);
        if currentTet_ids(3)<0
            currentTet_ids(3)=abs(currentTet_ids(3));
            Y3=Y.DataRow(currentTet_ids(3),:);
        else
            Y3 = Cell.FaceCentres.DataRow(currentTet_ids(3),:);
        end
        
        currentTet = [Y1; Y2; Y3; cellNuclei];
        %is currentTet correct??
    end
end
end

