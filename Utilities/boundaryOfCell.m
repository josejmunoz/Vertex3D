function [newVertOrder] = boundaryOfCell(verticesOfCell)
%BOUNDARYOFCELL Summary of this function goes here
% Here we consider 3 methods to connect the vertices, and we choose the
% method with more area into the polyshape.
% By Pedro J. Gomez Galvez, modified

    imaginaryCentroidMeanVert = mean(verticesOfCell);
    vectorForAngMean = bsxfun(@minus, verticesOfCell, imaginaryCentroidMeanVert );
    thMean = atan2(vectorForAngMean(:,2),vectorForAngMean(:,1));
    [~, angleOrderMean] = sort(thMean);
    %newVertOrderMean = [newVertOrderMean; newVertOrderMean(1,:)];

    %areaMeanCentroid = polyarea(newVertOrderMean(:,1),newVertOrderMean(:,2));
    
    newVertOrder = horzcat(angleOrderMean, vertcat(angleOrderMean(2:end), angleOrderMean(1)));

%     userConfig = struct('xy',verticesOfCell, 'showProg',false,'showResult',false);
%     resultStruct = tspo_ga(userConfig);
%     orderBoundary = [resultStruct.optRoute resultStruct.optRoute(1)];
%     
%     newVertSalesman = verticesOfCell(orderBoundary(1:end-1), :);
%     newVertSalesman = [newVertSalesman; newVertSalesman(1,:)];
%     areaVertSalesman = polyarea(newVertSalesman(:,1),newVertSalesman(:,2));
%     
%     if (areaVertSalesman >= areaMeanCentroid) && (areaVertSalesman >= areaCellCentroid)
%         newVertOrder = newVertSalesman;
%     else
%         if areaMeanCentroid >= areaCellCentroid
%             newVertOrder = newVertOrderMean;
%         else
%             newVertOrder = newVertOrderCent;
%         end
%     end
end

