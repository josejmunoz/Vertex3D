function [ quartetsOfNeighs ] = buildQuartetsOfNeighs2D(neighbours)

    quartetsOfNeighs={};
    for nCell=1:length(neighbours)
       
        neighCell=neighbours{nCell};
        interceptCells=cell(length(neighCell),1);
        for cellJ=1:length(neighCell)
            commonCells=intersect(neighCell,neighbours{neighCell(cellJ)});
            if length(commonCells)>2
               interceptCells{cellJ}=[commonCells' neighCell(cellJ) nCell]; 
            end
        end
        
        interceptCells=interceptCells(cellfun(@(x) ~isempty(x),interceptCells));
        
        if ~isempty(interceptCells)
            for indexA=1:length(interceptCells)-1
                for indexB=indexA+1:length(interceptCells)
                    intersectionCells = intersect(interceptCells{indexA},interceptCells{indexB});
                    if length(intersectionCells) >= 4
                        quartetsOfNeighs{end+1,1}=nchoosek(intersectionCells, 4);
                    end
                end
            end
        end
    end
    
    quartetsOfNeighs=unique(sort(cell2mat(quartetsOfNeighs), 2),'rows');
                
end