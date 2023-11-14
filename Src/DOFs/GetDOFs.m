function [Dofs]=GetDOFs(Geo, Set)
    % Define free and constrained vertices:
    %   1) Vertices with y-coordinates > Set.VPrescribed are those to be prescribed (pulled)
    %   2) Vertices with y-coordinates < Set.VFixed are those to be fixed
    %   3) the rest are set to be free
    % TODO FIXME HARDCODE
    dim = 3;
    gconstrained = zeros((Geo.numY+Geo.numF+Geo.nCells)*3, 1);
    gprescribed  = zeros((Geo.numY+Geo.numF+Geo.nCells)*3, 1);
    
    %% Fix border vertices
%     allIds = vertcat(Geo.Cells(Geo.BorderCells).globalIds);
%     allYs = vertcat(Geo.Cells(Geo.BorderCells).Y);
%     vertices_Top = allYs(:, 3) > Geo.CellHeightOriginal/2;
%     borderVertices_Top = allIds(vertices_Top);
%     vertices_Bottom = allYs(:, 3) < -Geo.CellHeightOriginal/2;
%     borderVertices_Bottom = allIds(vertices_Bottom);
    borderIds = vertcat(Geo.Cells(Geo.BorderCells).globalIds);
    %%
    for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        Y     = Geo.Cells(c).Y;
        gIDsY = Geo.Cells(c).globalIds;
        for f = 1:length(Geo.Cells(c).Faces)
            Face = Geo.Cells(c).Faces(f);
            if Face.Centre(2) < Set.VFixd
                gconstrained(dim*(Face.globalIds-1)+1:dim*Face.globalIds) = 1;
            elseif Face.Centre(2) > Set.VPrescribed
                gprescribed(dim*(Face.globalIds-1)+2) = 1;
				if Set.BC == 1
                	gconstrained(dim*(Face.globalIds-1)+1) = 1;
                	gconstrained(dim*(Face.globalIds-1)+3) = 1;
				end
            end
        end
        
%         if ismember(c, Geo.BorderCells)
%             fixY = ones(size(Y(:,2)));
%             preY = ones(size(Y(:,2)));
%         else 
            fixY = Y(:,2) < Set.VFixd | ismember(gIDsY, borderIds);
            preY = Y(:,2) > Set.VPrescribed | ismember(gIDsY, borderIds);
%         end
        for ff = 1:length(find(fixY))
            idx = find(fixY);
            idx = idx(ff);
            gconstrained(dim*(gIDsY(idx)-1)+1:dim*gIDsY(idx)) = 1;
        end

        for ff = 1:length(find(preY))
            idx = find(preY);
            idx = idx(ff);
            gprescribed(dim*(gIDsY(idx)-1)+1:dim*gIDsY(idx)) = 1;
        end
        
% 		if Set.BC == 1 % TODO FIXME Do not constrain this in compress...
%         	gconstrained(dim*(gIDsY(preY)-1)+1) = 1;
%         	gconstrained(dim*(gIDsY(preY)-1)+3) = 1;
% 		end
    end
    Dofs.Free = find(gconstrained==0 & gprescribed==0);
    Dofs.Fix  = [find(gconstrained); find(gprescribed)];
    Dofs.FixP = find(gprescribed);
    Dofs.FixC = find(gconstrained);
end