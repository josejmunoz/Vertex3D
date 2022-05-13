function [edges, edgesSharedbyCells] = BuildEdges(Tets, FaceIds, FaceCentre, X, Ys, nonDeadCells)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BuildEdges:										  
    %   Obtain the local ids of the edges that define a Face. The order of 
    %   such indices is reordered to produce positive areas and volumes.
    % Input:															  
    %   Tets		: All tetrahedras
    %   FaceIds		: Local indices defining the tetrahedras in the face
    %   FaceCentre	: Centre of the Face
    %   X           : Centre of the Cell containing the Face
    %   Ys          : All Vertices of the Cell
    % Output:															  
    %   edges       : Local indices of the vertices forming the 
    %	face. That is Geo.Cells(c).Y(edges(e,:),:) will give vertices
    %   defining the edge. Used also for triangle computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	FaceTets = Tets(FaceIds,:);
	tet_order = zeros(length(FaceTets),1);
	% TODO FIXME, initialize, there was a bug here. Is there a more
	% clean way to write it ?
	tet_order(1) = 1;
	prev_tet  = FaceTets(1,:);
	%% Tetrahedra ordering
	% Order the tetrahedras defining the vertices of the face in a general
	% clock-wise manner. 
    if size(FaceTets,1) > 3
        for yi = 2:length(FaceTets)
		    i = sum(ismember(FaceTets, prev_tet),2)==3;
		    i = i & ~ismember(1:length(FaceTets),tet_order)';
		    i = find(i);
            if isempty(i)
                edges = [];
                return
            end
		    tet_order(yi) = i(1);
		    prev_tet = FaceTets(i(1),:);
        end
    else
        % TODO FIXME is this enough??? will it get flipped later if not
        % correct ???
        tet_order = [1 2 3]';
    end
	surf_ids  = 1:length(Tets); 
	surf_ids  = surf_ids(FaceIds);
    if length(surf_ids) < 3
		% Something went really wrong in the simulation, or a flip being 
		% tested would result in another face being just an edge.
        edges = [];
        return
    end
	surf_ids  = surf_ids(tet_order);
    % TODO FIXME IS THIS ACCEPTABLE?
%     if size(FaceTets,1) > 3
	%% Vertices ordering
	% The clockwise ordering might be incorrect for some cases, which need
	% reordering.
	if size(FaceTets,1) == 3
		centre = sum(Ys(surf_ids,:))/3;
	else
		centre = FaceCentre;
	end
    Order=0;
    for iii=1:length(surf_ids)
	    if iii==length(surf_ids)
		    v1=Ys(surf_ids(iii),:)-centre;
		    v2=Ys(surf_ids(1),:)-centre;
		    Order=Order+dot(cross(v1,v2),centre-X)/length(surf_ids);
	    else
		    v1=Ys(surf_ids(iii),:)-centre;
		    v2=Ys(surf_ids(iii+1),:)-centre;
		    Order=Order+dot(cross(v1,v2),centre-X)/length(surf_ids);
	    end
    end
	if Order<0 
	    surf_ids=flip(surf_ids);
    end
	%% Build edges and identify the ones shared by different cells
	edges = zeros(length(surf_ids), 2);
    edgesSharedbyCells = zeros(length(surf_ids), 1);
	for yf = 1:length(surf_ids)-1
		edges(yf,:) = [surf_ids(yf) surf_ids(yf+1)];
        %edges shared by different cells
        edgesSharedbyCells(yf) = sum(ismember(Tets(edges(yf, 1) , :), nonDeadCells)) >= 2 & sum(ismember(Tets(edges(yf, 2) , :), nonDeadCells)) >= 2;
	end
	edges(end,:) = [surf_ids(end) surf_ids(1)];
    edgesSharedbyCells(end) = sum(ismember(Tets(edges(end, 1) , :), nonDeadCells)) >= 2 & sum(ismember(Tets(edges(end, 2) , :), nonDeadCells)) >= 2;
    
end