function [Twg_upsampled, X, layerNodeIDs] = upsampleTetMesh(Twg, X, layerNodeIDs)
%UPSAMPLETETMESH Summary of this function goes here
%   Based on: https://github.com/alecjacobson/gptoolbox/blob/master/mesh/upsample.m
% 
%   %      o           o     
%   %     / \         / \     
%   %    x   x  ---> o---o    
%   %   /     \     / \ / \  
%   %  o---x---o   o---o---o  
%   function [U14,F14,E14] = one_four(offset,V,F)
%     % compute midpoints (actually repeats, one midpoint per edge per face)
%     E14 = [F(:,2) F(:,3);F(:,3) F(:,1);F(:,1) F(:,2)];
%     U14 = (V(E14(:,1),:)+V(E14(:,2),:))/2;
%     % indices of midpoints
%     nu = size(U14,1);
%     i1 = offset+(1:(nu/3))';
%     i2 = offset+((nu/3)) + (1:(nu/3))';
%     i3 = offset+((nu/3)+(nu/3)) + (1:(nu/3))';
%     % new face indices, 4 new faces for each original face. As if we simply
%     % ignored the duplicates in m and had appended m to V
%     F14 = [ F(:,1) i3 i2 ; F(:,2) i1 i3 ; F(:,3) i2 i1 ; i1 i2 i3];
%   end


tetsToUpsample = Twg(sum(ismember(Twg, layerNodeIDs), 2) >= 3, :);

Twg_upsampled = zeros(size(tetsToUpsample, 1)*3, 4);

lastID = size(X, 1);
for nTet = 1:size(tetsToUpsample, 1)
    tet = tetsToUpsample(nTet, :);
    gNodes = tet(ismember(tet, layerNodeIDs));
    mainNodes = setdiff(tet, gNodes);

    newNode = mean(X(gNodes, :));
    newNodeID = lastID + 1;
    lastID = newNodeID;

    X(newNodeID, :) = newNode;
    layerNodeIDs(end+1) = newNodeID;


    newTriangles = nchoosek([gNodes newNodeID], 3);
    
    newTriangles(all(ismember(newTriangles, gNodes), 2), :) = [];
    Twg_upsampled(3*(nTet-1)+1:3*(nTet), :) = horzcat(newTriangles, ones(3, 1) * mainNodes);
end

Twg(sum(ismember(Twg, layerNodeIDs), 2) >= 3, :) = [];
new_Twg_upsampled = vertcat(Twg, Twg_upsampled);

oldLayerNodeIDs = layerNodeIDs;
numberOfNewTets = 2;
tetsToUpsample = new_Twg_upsampled(sum(ismember(new_Twg_upsampled, layerNodeIDs), 2) == numberOfNewTets, :);
Twg_upsampled = zeros(size(tetsToUpsample, 1)*numberOfNewTets, 4);
Twg_upsampled_Adapted = [];

tetsToRemove = [];

for nTet = 1:size(tetsToUpsample, 1)
    tet = tetsToUpsample(nTet, :);
    gNodes = tet(ismember(tet, layerNodeIDs));
    mainNodes = setdiff(tet, gNodes);

    newNode = mean(X(gNodes, :));
    newNodeID = lastID + 1;
    lastID = newNodeID;

    X(newNodeID, :) = newNode;
    layerNodeIDs(end+1) = newNodeID;
    
    newTriangles = nchoosek([gNodes newNodeID], numberOfNewTets);
    newTriangles(all(ismember(newTriangles, gNodes), 2), :) = [];
    
    Twg_upsampled(numberOfNewTets*(nTet-1)+1:numberOfNewTets*(nTet), :) = horzcat(newTriangles, repmat(mainNodes, [numberOfNewTets, 1]));

    tetsToAdaptDueToChange = new_Twg_upsampled(sum(ismember(new_Twg_upsampled, gNodes), 2) >= numberOfNewTets, :);
    tetsToRemove = vertcat(tetsToRemove, tetsToAdaptDueToChange);
    tetsToAdaptDueToChange(all(ismember(tetsToAdaptDueToChange, tet), 2), :) = [];
    
    for nTetToAdapt = 1:size(tetsToAdaptDueToChange, 1)
        tetToAdapt = tetsToAdaptDueToChange(nTetToAdapt, :);

        nodesToConnect = tetToAdapt(~ismember(tetToAdapt, gNodes));

        Twg_upsampled_Adapted(end+1:end+2, :) = horzcat(newTriangles, repmat(nodesToConnect, [numberOfNewTets, 1]));
    end
end

new_Twg_upsampled(ismember(new_Twg_upsampled, tetsToRemove, 'rows'), :) = [];
new_Twg_upsampled = vertcat(new_Twg_upsampled, Twg_upsampled);
new_Twg_upsampled = vertcat(new_Twg_upsampled, Twg_upsampled_Adapted);

Twg_upsampled = new_Twg_upsampled;

end

