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
    allIds = vertcat(Geo.Cells(1:Geo.nCells).globalIds);
    allYs = vertcat(Geo.Cells(1:Geo.nCells).Y);
    vertices_Top = allYs(:, 3) > Set.CellHeight/4 & allYs(:, 3) < Set.CellHeight;
    vertices_TopIds = allIds(vertices_Top);
    borderVertices_Top = vertices_TopIds(boundary(allYs(vertices_Top, 1:2), 0.7));
    vertices_Bottom = allYs(:, 3) < -Set.CellHeight/4 & allYs(:, 3) > -Set.CellHeight;
    vertices_BottomIds = allIds(vertices_Bottom);
    borderVertices_Bottom = vertices_BottomIds(boundary(allYs(vertices_Bottom, 1:2), 0.7));
    %%
    for c = 1:Geo.nCells
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
        fixY = Y(:,2) < Set.VFixd | ismember(gIDsY, borderVertices_Top) | ismember(gIDsY, borderVertices_Bottom);
        preY = Y(:,2) > Set.VPrescribed | ismember(gIDsY, borderVertices_Top) | ismember(gIDsY, borderVertices_Bottom);
        for ff = 1:length(find(fixY))
            idx = find(fixY);
            idx = idx(ff);
            gconstrained(dim*(gIDsY(idx)-1)+1:dim*gIDsY(idx)) = 1;
        end
        
        gprescribed(dim*(gIDsY(preY)-1)+1:dim*(gIDsY(preY)-1)) = 1;
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