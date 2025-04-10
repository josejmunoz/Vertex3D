function [Tris] = BuildEdges(Tets, FaceIds, FaceCentre, FaceInterfaceType, X, Ys, nonDeadCells)
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
	TrisFields = ["Edge", "Area", "AspectRatio", "EdgeLength", "LengthsToCentre", "SharedByCells", "Location", "ContractileG", "ContractilityValue", "EdgeLength_time", "pastContractilityValue"]; 
    Tris = BuildStructArray(sum(FaceIds), TrisFields);

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
                valid_orders = find_valid_tetrahedra_orders(Tets, FaceIds);
                if isempty(valid_orders) == false
                    tet_order = valid_orders{1};
                    break
                else
                    ME = MException('BuildEdges:TetrahedraOrdering', ... 
                        sprintf('Cannot create a face with these tetrahedra'));
                    throw(ME);
                end
            end
	        tet_order(yi) = i(1);
	        prev_tet = FaceTets(i(1),:);
        end
        % Last one should match with the first one
        if sum(ismember(FaceTets(1, :), prev_tet),2)~=3
            valid_orders = find_valid_tetrahedra_orders(Tets, FaceIds);
            if isempty(valid_orders) == false
                tet_order = valid_orders{1};
            else
                ME = MException('BuildEdges:TetrahedraOrdering', ...
                    sprintf('Cannot create a face with these tetrahedra'));
                throw(ME);
            end
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
        ME = MException('BuildEdges:TetrahedraMinSize', ...
            'Length of the face is lower than 3');
        throw(ME);
    end
	surf_ids  = surf_ids(tet_order);
	%% Vertices ordering
	% The clockwise ordering might be incorrect for some cases, which need
	% reordering.
    Order=zeros(1, length(surf_ids));
    for iii=1:length(surf_ids)
	    if iii==length(surf_ids)
		    v1=Ys(surf_ids(iii),:)-FaceCentre;
		    v2=Ys(surf_ids(1),:)-FaceCentre;
		    Order(iii)=dot(cross(v1,v2),FaceCentre-X)/length(surf_ids);
	    else
		    v1=Ys(surf_ids(iii),:)-FaceCentre;
		    v2=Ys(surf_ids(iii+1),:)-FaceCentre;
		    Order(iii)=dot(cross(v1,v2),FaceCentre-X)/length(surf_ids);
	    end
    end
	if all(Order<0) 
	    surf_ids=flip(surf_ids);
%     elseif any(Order<0)
%         disp('possible error')
%         surf_ids(Order < 0)=flip(surf_ids(Order < 0));
    end
    
	%% Build edges and identify the ones shared by different cells
	for currentTri = 1:length(surf_ids)-1
		Tris(currentTri).Edge = [surf_ids(currentTri) surf_ids(currentTri+1)];
        %edges shared by different cells
        currentTris_1 = Tets(Tris(currentTri).Edge(1), :);
        currentTris_2 = Tets(Tris(currentTri).Edge(2), :);
        Tris(currentTri).SharedByCells = intersect(currentTris_1(ismember(currentTris_1, nonDeadCells)), currentTris_2(ismember(currentTris_2, nonDeadCells)));
        
        % Compute Tris aspect ratio, edge length and LengthsToCentre
        [Tris(currentTri).EdgeLength, Tris(currentTri).LengthsToCentre, Tris(currentTri).AspectRatio] = ComputeTriLengthMeasurements(Tris, Ys, currentTri, FaceCentre);
        Tris(currentTri).EdgeLength_time = [0, Tris(currentTri).EdgeLength];
	end
	Tris(length(surf_ids)).Edge = [surf_ids(end) surf_ids(1)];
    currentTris_1 = Tets(Tris(length(surf_ids)).Edge(1), :);
    currentTris_2 = Tets(Tris(length(surf_ids)).Edge(2), :);
    Tris(length(surf_ids)).SharedByCells = intersect(currentTris_1(ismember(currentTris_1, nonDeadCells)), currentTris_2(ismember(currentTris_2, nonDeadCells)));
    
    % Compute Tris aspect ratio, edge length and LengthsToCentre
    [Tris(length(surf_ids)).EdgeLength, Tris(length(surf_ids)).LengthsToCentre, Tris(length(surf_ids)).AspectRatio] = ComputeTriLengthMeasurements(Tris, Ys, length(surf_ids), FaceCentre);
    Tris(length(surf_ids)).EdgeLength_time = [0, Tris(length(surf_ids)).EdgeLength];

    % Compute Tris area
    [~, triAreas] = ComputeFaceArea(vertcat(Tris.Edge), Ys, FaceCentre);
    [Tris.Area] = triAreas{:};
    
    % Add edge location: 'Top/Bottom/Lateral'
    [Tris.Location] = deal(FaceInterfaceType);
    
    % Initialize forces
    [Tris.ContractileG] = deal(0);
end

function valid_orders = find_valid_tetrahedra_orders(Tets, FaceIds) 
    FaceTets = Tets(FaceIds, :); 
    n_tets = size(FaceTets, 1); 
    valid_orders = {};  
    
    for iter = 1:100000000
        perm = randperm(n_tets);
        is_valid = true;
        for i = 2:n_tets
            prev_tet = FaceTets(perm(i-1), :);
            curr_tet = FaceTets(perm(i), :);
            if sum(ismember(curr_tet, prev_tet)) ~= 3
                is_valid = false;
                break;
            end
        end
        if is_valid && sum(ismember(FaceTets(perm(1), :), FaceTets(perm(end), :))) == 3
            valid_orders{end+1} = perm;
            break
        end
    end
end