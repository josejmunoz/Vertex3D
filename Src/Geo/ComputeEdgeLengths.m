function [edgeLengths] = ComputeEdgeLengths(Cell)
%% Compute the length of the segments between vertices Xs
% loop on cells
allEdges = [];
for numCell=1:obj.n
    obj.EdgeLengths{numCell}=zeros(size(obj.Cv{numCell},1),1);
    % loop on edges
    for e=1:size(obj.Cv{numCell},1)
        %                Cv(1,:) ---- >always vertex
        Y1=Y.DataRow(obj.Cv{numCell}(e,1),:);
        if  obj.Cv{numCell}(e,2) > 0
            %                    Cv(2,:)>0 ---- > is vertex
            Y2=Y.DataRow(obj.Cv{numCell}(e,2),:);
        else
            %                    Cv(2,:)<0 ---- > is face center
            Y2=obj.FaceCentres.DataRow(abs(obj.Cv{numCell}(e,2)),:);
        end
        % Compute Length
        allEdges = vertcat(allEdges, sort([obj.Cv{numCell}(e,1) obj.Cv{numCell}(e,2)]));
        obj.EdgeLengths{numCell}(e)=norm(Y1-Y2);
    end
    
    obj.ContractileForces{numCell} = zeros(size(obj.Cv{numCell}, 1), 1);
end
[~, uniqueEdges] = unique(allEdges, 'rows');
end