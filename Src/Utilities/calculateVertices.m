function [ verticesInfo ] = calculateVertices( labelledImg, neighbours, ratio)
%CALCULATEVERTICES Summary of this function goes here
%   With a labelled image as input, the objective is get all vertex for
%   each cell
%   Developed by Pedro Gomez Galvez

    se=strel('disk',ratio);
    neighboursVertices = buildTripletsOfNeighs( neighbours );%intersect dilatation of each cell of triplet
    vertices = cell(size(neighboursVertices, 1), 1);
        
    % We first calculate the perimeter of the cell to improve efficiency
    % If the image is small, is better not to use bwperim
    % For larger images it improves a lot the efficiency
    dilatedCells=cell(max(max(labelledImg)),1);
    
    for i=1:max(max(labelledImg))
        BW=zeros(size(labelledImg));
        BW(labelledImg==i)=1;
        BW_dilated=imdilate(bwperim(BW),se);
        dilatedCells{i}=BW_dilated;
    end
    
    %the overlapping between labelledImg cells will be the vertex
    borderImg=zeros(size(labelledImg));
    borderImg(labelledImg>-1)=1;
    for numTriplet = 1 : size(neighboursVertices,1)
              
        BW1_dilate=dilatedCells{neighboursVertices(numTriplet, 1),1};
        BW2_dilate=dilatedCells{neighboursVertices(numTriplet, 2),1};
        BW3_dilate=dilatedCells{neighboursVertices(numTriplet, 3),1};
        
        %It is better use '&' than '.*' in this function
        [row,col]=find((BW1_dilate.*BW2_dilate.*BW3_dilate.*borderImg)==1);
        
        if length(row)>1
            if ~ismember(round(mean(col)),col)
                vertices{numTriplet,1}=round(mean([row(col > mean(col)),col(col > mean(col))]));
                vertices{numTriplet,2}=round(mean([row(col < mean(col)),col(col < mean(col))]));
            else
                vertices{numTriplet} = round(mean([row,col]));
            end
        else
            vertices{numTriplet} = [row,col];
        end
        
    end
    
    %storing vertices and deleting artefacts
    verticesInfo.location = vertices;
    verticesInfo.connectedCells = neighboursVertices;
    
    notEmptyCells=cellfun(@(x) ~isempty(x),verticesInfo.location,'UniformOutput',true);
    if size(verticesInfo.location,2)==2
        verticesInfo.location=[verticesInfo.location(notEmptyCells(:,1),1);verticesInfo.location(notEmptyCells(:,2),2)];
        verticesInfo.connectedCells=[verticesInfo.connectedCells(notEmptyCells(:,1),:);verticesInfo.connectedCells(notEmptyCells(:,2),:)];
    else
        verticesInfo.location=verticesInfo.location(notEmptyCells,:);
        verticesInfo.connectedCells=verticesInfo.connectedCells(notEmptyCells,:);
    end
    
    verticesInfo.location = vertcat(verticesInfo.location{:});

end