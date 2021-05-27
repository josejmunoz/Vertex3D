classdef CellClass
    %% Cell class
    properties
        Int           % - Cell nodes (array-structure ,  Size={1 NumCells}):
        %       IDs of internal nodes (Cell centres)
        %---------------------------------------------------------------------
        Tris                  % -Cell-surface Triangulation (Type=cell-structure ,  Size={NumCells 1})
        %     Each cell is an array of size [ntris 3] with triangles defining cell surfaces.
        %     e.g. Cell.Tris{1}(2,:)=[v1 v2 v3] -> such that v1, v2 and v3 are the vertices which
        %     defines Triangle (2) of  Cell (1).
        %     Remark:- v1 and v2 always correspond to vertices (their position can be found as Y.DataRow([v1 v2],:))
        %            - v3 >0 correspond to a face-centre  (its position can be found as Cell.FaceCentres.DataRow(v3,:))
        %            - v3 <0 correspond to a vertex (its position can be found as Y.DataRow(abs(v3),:))
        %--------------------------------------------------------------------
        cTet                  % -Connected tetrahedrons (Type=cell-structure ,  Size={NumCells 1}):
        %     Each cell (Cell.cTet{i}) is an array of size [ntet 4] with nodal tetrahedrons connected to the node of Cell i.
        %--------------------------------------------------------------------
        cTetID                %% -Connected tetrahedrons IDs (Type=cell-structure ,  Size={NumCells 1}):
        %     Each cell (Cell.cTetID{i}) is an array of size [1 ntet] with IDs of nodal tetrahedrons connected to the node of Cell i.
        %    (IDs-> their index in T).
        %--------------------------------------------------------------------
        cNodes                % -Connected Nodes (Type=cell-structure ,  Size={NumCells 1}):
        %  Each cell (Cell.cNodes{i})is an array of size [1 nNodes] with IDs of nodes connected to the node of Cell i.
        %--------------------------------------------------------------------
        Cv                    % -Edges (connectivity of vertices) (Type=cell-structure ,  Size={NumCells 1}):
        %    Each cell (Cell.Cv{i})is an array of size [nEdges 2] with all the edges between each two vertices of Cell i .
        %     e.g. Cell.Cv{1}(2,:)=[v1 v2] -> such that v1 and v2 are the vertices that correspond to edge (2) of  Cell (1).
        %     Remark:- v1 always correspond to a vertex   (its position can be found as Y.DataRow(v1,:))
        %            - v2 >0  also correspond to a vertex (its position can be found as Y.DataRow(v2,:))
        %            - v2 <0  correspond to a face-centre  (its position can be found as Cell.FaceCentres.DataRow(abs(v2),:))
        %--------------------------------------------------------------------
        CvID                 % -The IDs of edges (connectivity of vertices) (Type=cell-structure ,  Size={NumCells 1}):
        %                        Each cell (Cell.CvID{i})is an array of size [nEdges 1] with the IDs of edges between each two vertices of Cell i
        %--------------------------------------------------------------------
        EdgeLengths              % -The length of Edges at the previous time-step (Type=cell-structure ,  Size={NumCells 1}):
        %                                 Each cell (Cell.EdgeLengths{i})is an array of size [nEdges 1] with the length of the edges between each two vertices.
        %--------------------------------------------------------------------
        EdgeLengthsn              % -The length of Edges at the previous time-step  (Type=cell-structure ,  Size={NumCells 1}):
        %                            Each cell (Cell.EdgeLengthn{i})is an array of size [nEdges 1] with the length of the edges between each two vertices.
        %--------------------------------------------------------------------
        EdgeLengths0_average         % -The Initial\reference length of Edges at the  (Type=cell-structure ,  Size={NumCells 1}):
        %                            Each cell (Cell.EdgeLengths0{i})is an array of size [nEdges 1] with the length of the edges between each two vertices.
        
        %--------------------------------------------------------------------
        FaceCentres          %% -Face centres  (Type=array-structure ,  Size={1 NumFaces}):
        % the x-y-z coordinates of face centres
        %
        FaceCentres0
        %--------------------------------------------------------------------
        Centre
        Centre0
        %--------------------------------------------------------------------
        nTotalTris           %% - Total number of triangles (Type=scalar)
        %--------------------------------------------------------------------
        n                    %% - Total number of Cells (Type=scalar)
        %--------------------------------------------------------------------
        Vol                  %% - Volume of Cells                   (Type=array-structure ,  Size=[NumCells1 ]):
        %--------------------------------------------------------------------
        Vol0                 %% - Initial\reference volume of Cells (Type=array-structure ,  Size=[NumCells1 ]):
        %--------------------------------------------------------------------
        SArea                %% - Surface area of Cells             (Type=array-structure ,  Size=[NumCells1 ]):
        %--------------------------------------------------------------------
        SArea0               %% - Initial\reference area of Cells   (Type=array-structure ,  Size=[NumCells1 ]):
        %--------------------------------------------------------------------
        Faces                %% - Cell Faces (Type=cell-structure ,  Size={NumCells 1}):
        %     - Cell.Faces{i}.nFaces          -> number of faces of Cell i            (Type=scalar)
        %     - Cell.Faces{i}.FaceCentresID   -> The IDs of faces of Cell i           (Type=array-structure, Size=[Cell{i}.nFaces 1])
        %     - Cell.Faces{i}.Vertices        -> The vertices of each face of Cell i  (Type=cell-structure, Size={Cell{i}.nFaces 1})
        %           (e.g. Cell.Faces{i}.Vertices{j} --> is an array the list the vertices of Face j of Cell i )
        %     - Cell.Faces{i}.Tris            -> The vertex-triangles of each face of Cell i  (Type=cell-structure, Size={Cell{i}.nFaces 1})
        %           (e.g. Cell.Faces{i}.Tris{j}(k,:)=[v1 v2 v3] --> v1, v2 and v3 are the vertices of triangle k of Face j of Fell k)
        %            Remark:- v1 and v2 always correspond to vertices (their position can be found as Y.DataRow([v1 v2],:))
        %                   - v3 >0 correspond to a face-centre  (its position can be found as Cell.Face Centres .DataRow(v3,:))
        %                   - v3 <0 correspond to a vertex (its position can be found as Y.DataRow(abs(v3),:))
        %--------------------------------------------------------------------
        SAreaTri              % - The area of Surfaces triangles  (Type=cell-structure ,  Size={NumCells 1}):
        %        Each cell (Cell.SAreaTri{i})is an array with the area of the triangles of Cell i.
        %--------------------------------------------------------------------
        SAreaTrin             % - The area of Surfaces triangles at the previous time-step (Type=cell-structure ,  Size={NumCells 1}):
        %        Each cell (Cell.SAreaTrin{i})is an array with the area of the triangles of Cell i at the previous time-step.
        %--------------------------------------------------------------------
        SAreaFace             % - The area of Faces  (Type=cell-structure ,  Size={NumCells 1}):
        %        Each cell (Cell.SAreaFace{i})is an array with the area of the all faces of Cell i.
        %--------------------------------------------------------------------
        SAreaFace0            % - The Initial\reference area of Faces  (Type=cell-structure ,  Size={NumCells 1}):
        %        Each cell (Cell.SAreaFaceo{i})is an array with the initial\reference area of the all faces of Cell i.
        %--------------------------------------------------------------------
        AssembleAll           % - Assembly condition  (logical)
        %       - Cell.AssembleAll = true  --> The contribution of all Cells is assembled in K and g
        %       - Cell.AssembleAll = false --> Only the contribution of Cell.Assemble is to be assembled in K and g
        %--------------------------------------------------------------------
        AssembleNodes         % - The nodes\cell centres of cells to be assembles in K and g. (Type=cell-structure)
        %--------------------------------------------------------------------
        RemodelledVertices    % - The new vertices which have just been remodelled (those to be assembled in the local problem)
        %--------------------------------------------------------------------
        Edges                 % Edges (Type=cell-structure ,  Size={NumCells 1}):
        %    Each cell (Cell.Edges{i})is an array of size [nEdges 4] with each row corresponds to an edge  all the edges between each two vertices.
        %    (e.g. Cell.Edges{i}(j,:)=[v1 v2 v3 v4] --> v1, v2, v3 and v4 are the two vertices of edge
        %      (j) besides the two vertices of the triangles that share edge (j))
        %     Remark:- v1,v2,v3 and v4 <=Y.n  correspond to a vertex         (its position can be found as Y.DataRow(v1,:))
        %            - v1,v2,v3 and v4 >Y.n   correspond to a face-centre    (its position can be found as Cell.FaceCentres.DataRow(v1-T.n,:))
        %--------------------------------------------------------------------
        EdgeLocation
        %--------------------------------------------------------------------
        ApicalVertices
        BasalVertices
        %--------------------------------------------------------------------
        ApicalBorderVertices
        BasalBorderVertices
        %
        BorderVertices % of the tissue
        %--------------------------------------------------------------------
        SubstrateForce
        %--------------------------------------------------------------------
        DebrisCells            % Debris cells
        %--------------------------------------------------------------------
        BorderCells
        %--------------------------------------------------------------------
        ContractileForces
    end
    methods
        function Cell = CellClass(nC,xInternal)
            if nargin > 0
                Cell.Int=xInternal;
                Cell.Tris=cell(nC,1);
                Cell.cTet=cell(nC,1);
                Cell.cTetID=cell(nC,1);
                Cell.cNodes=cell(nC,1);
                Cell.Cv=cell(nC,1);
                Cell.CvID=cell(nC,1);
                Cell.n=length(xInternal);
                Cell.Vol=zeros(Cell.n,1);
                Cell.Vol0=zeros(Cell.n,1);
                Cell.SArea=zeros(Cell.n,1);
                Cell.SArea0=zeros(Cell.n,1);
                Cell.SAreaTri=cell(Cell.n,1);
                Cell.SAreaTrin=cell(Cell.n,1);
                Cell.SAreaFace=cell(Cell.n,1);
                Cell.SAreaFace0=cell(Cell.n,1);
                Cell.Faces=cell(nC,1);
                Cell.AssembleAll=true;
                Cell.AssembleNodes=[];
                Cell.RemodelledVertices=[];
                Cell.Edges=cell(nC,1);
                Cell.EdgeLengths=cell(nC,1);
                Cell.EdgeLengths0_average=-1;
                Cell.EdgeLengthsn=cell(nC,1);
                Cell.DebrisCells=false(nC, 1);
                Cell.ContractileForces=cell(nC, 1);
                Cell.EdgeLocation=cell(nC,1);
                Cell.ApicalVertices=cell(nC, 1);
                Cell.BasalVertices=cell(nC, 1);
                Cell.SubstrateForce=cell(nC, 1);
                Cell.ApicalBorderVertices = cell(nC, 1);
                Cell.BasalBorderVertices = cell(nC, 1);
                Cell.BorderVertices = [];
                Cell.BorderCells=false(nC, 1);
                Cell.Centre = zeros(nC, 3);
                Cell.Centre0 = zeros(nC, 3);
            end
        end
        
        %% Ablate cells
        % We assume there will only be one ablation. Thus, we could remove
        % the IDs from the middle.
        % TODO: Check if Cell.Int need to be in order in consecutive numbers
        function Cell = AblateCells(obj, cellsToRemove)
            obj.DebrisCells(ismember(obj.Int, cellsToRemove)) = true;
            Cell = obj;
        end 
        
        function Cell = removeCells(obj, cellsToRemove)
            obj.Int(cellsToRemove) = [];
            
            %Consequences:
            %cell structures
            obj.Tris(cellsToRemove) = [];
            obj.cTet(cellsToRemove) = [];
            obj.cTetID(cellsToRemove) = [];
            obj.cNodes(cellsToRemove) = [];
            obj.Cv(cellsToRemove) = [];
            obj.CvID(cellsToRemove) = [];
            
            %Update info of the cells
            obj.Vol(cellsToRemove) = [];
            obj.Vol0(cellsToRemove) = [];
            obj.SArea(cellsToRemove) = [];
            obj.SArea0(cellsToRemove) = [];
            
            % Update info of the cells, with cell structure
            obj.SAreaTri(cellsToRemove) = [];
            obj.SAreaTrin(cellsToRemove) = [];
            obj.SAreaFace(cellsToRemove) = [];
            obj.SAreaFace0(cellsToRemove) = [];
            obj.Faces(cellsToRemove) = [];
            obj.Edges(cellsToRemove) = [];
            obj.EdgeLengths(cellsToRemove) = [];
            obj.EdgeLengthsn(cellsToRemove) = [];
            obj.ContractileForces(cellsToRemove) = [];
            obj.DebrisCells(cellsToRemove) = [];
            
            obj.n = obj.n - sum(cellsToRemove);
            obj.nTotalTris = sum(cellfun(@(X) size(X,1), obj.Tris));
            
            Cell = obj;
        end
        
        %% Compute the length of the segments between vertices Xs
        function [obj, uniqueEdges] = computeEdgeLengths(obj, Y)
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
        
        %%
        function [obj, featuresTable, resultingImage] = exportTableWithCellFeatures(obj, Y, timeStep, Faces, Set)
            %% Features to obtain per tissue:
            % Avg Cell height
            
            
            %% Features to obtain per cell:
            cellCellFaces = find(Faces.InterfaceType == 1);
            
            featuresTable_cell = [];
            for numCell = obj.Int
                % - Cell height
                % avg and std distance between connected apical and basal
                % vertices
                cellHeight = mean(obj.EdgeLengths{numCell}(obj.EdgeLocation{numCell} == 1));
                cellHeightSTD = std(obj.EdgeLengths{numCell}(obj.EdgeLocation{numCell} == 1));
                % - Volume
                cellVolume = obj.Vol(numCell);

                % - Apical area,  basal area and lateral area
                cellSurfaceArea = obj.SArea(numCell);

                trianglesArea = obj.SAreaTri{numCell};
                apicalFaceCentres = abs(obj.ApicalVertices{numCell}(obj.ApicalVertices{numCell} < 0));
                apicalTriangles = all(ismember(obj.Tris{numCell}(:, 1:2), obj.ApicalVertices{numCell}), 2) & ismember(obj.Tris{numCell}(:, 3), apicalFaceCentres);
                apicalArea = sum(trianglesArea(apicalTriangles));
                basalFaceCentres = abs(obj.BasalVertices{numCell}(obj.BasalVertices{numCell} < 0));
                basalTriangles = all(ismember(obj.Tris{numCell}(:, 1:2), obj.BasalVertices{numCell}), 2) & ismember(obj.Tris{numCell}(:, 3), basalFaceCentres);
                basalArea = sum(trianglesArea(basalTriangles));
                
                lateralTrianglesArea = cellSurfaceArea - apicalArea - basalArea;
                
                % - Shared lateral area per neighbour (avg and std)
                currentFaceIDs = obj.Faces{numCell}.FaceCentresID;
                cellCellAreas = obj.SAreaFace{1}(ismember(currentFaceIDs, cellCellFaces));
                lateralAreaSharedPerNeighbour_AVG = mean(cellCellAreas);
                lateralAreaSharedPerNeighbour_STD = std(cellCellAreas);
                
                % - Neighbours: 3D neighbours, apical and basal neighbours,
                % polygon distribution
                apicalNeighbours = [];
                basalNeighbours = [];
                
                basalBorderVerticesIDs = obj.BasalVertices{numCell}(obj.BasalBorderVertices{numCell});
                apicalBorderVerticesIDs = obj.ApicalVertices{numCell}(obj.ApicalBorderVertices{numCell});
                for neighbourCell = 1:obj.n
                    if neighbourCell ~= numCell
                        %Apical neighbours
                        if any(ismember(apicalBorderVerticesIDs, obj.ApicalVertices{neighbourCell}(obj.ApicalBorderVertices{neighbourCell})))
                            apicalNeighbours = [apicalNeighbours, neighbourCell];
                        end
                        
                        %Basal neighbours
                        if any(ismember(basalBorderVerticesIDs, obj.BasalVertices{neighbourCell}(obj.BasalBorderVertices{neighbourCell})))
                            basalNeighbours = [basalNeighbours, neighbourCell];
                        end
                    end
                end
                
                % - Convexity/concavity of cell
                % - info normalized: 35microns cell height initial
                Set.CellHeight;
                
                featuresTable_cell = [featuresTable_cell; numCell, obj.BorderCells(numCell), obj.DebrisCells(numCell), cellHeight, ...
                    cellHeightSTD, cellVolume, cellSurfaceArea, apicalArea, ...
                    basalArea, lateralTrianglesArea, lateralAreaSharedPerNeighbour_AVG, ...
                    lateralAreaSharedPerNeighbour_STD, {apicalNeighbours}, ...
                    numel(apicalNeighbours), {basalNeighbours}, numel(basalNeighbours)];
            end
            
            featuresTable = cell2table(featuresTable_cell, 'VariableNames', ...
                {'ID', 'IsBorderCell', 'IsDebrisCell', 'HeightAVG', 'HeightSTD', ...
                'Volume', 'SurfaceArea', 'ApicalArea', 'BasalArea', ...
                'LateralArea', 'LateralAreaSharedPerNeighbourAVG', ...
                'LateralAreaSharedPerNeighbourASTD', 'ApicalNeighbours', ...
                'NumApicalNeighbours', 'BasalNeigbhours', 'NumBasalNeighbours'});
            
        end
        
        %%
        function [obj] = computeEdgeLocation(obj, Y)
            for numCell=1:obj.n
                currentEdgesOfCell = obj.Cv{numCell};
                uniqueCurrentVertices = unique(currentEdgesOfCell(currentEdgesOfCell > 0));
                if obj.DebrisCells(numCell)
                    remainingEdges = vertcat(obj.Cv{setdiff(find(obj.DebrisCells == 0), numCell)});
                else
                    remainingEdges = vertcat(obj.Cv{setdiff(1:obj.n, numCell)});
                end
                uniqueRemainingEdges = unique(remainingEdges(remainingEdges > 0));

                sharedVertices = uniqueCurrentVertices(ismember(uniqueCurrentVertices, uniqueRemainingEdges));
                
                edgesToFaces = ismember(currentEdgesOfCell(:, 1), sharedVertices) & currentEdgesOfCell(:, 2) < 0;
                [numElements, elements] = hist(currentEdgesOfCell(edgesToFaces, 2), unique(currentEdgesOfCell(edgesToFaces, 2)));
                cellCellFaceCentres = elements(numElements > 3);
                midZ = obj.FaceCentres.DataRow(abs(cellCellFaceCentres),3);
                
                if length(midZ) <= 3
                    %It may be a border cell
                end
                
                apicoBasalVertices = Y.DataRow(sharedVertices,:);
                upperVerticesBorder = sharedVertices(apicoBasalVertices(:, 3) > mean(midZ));
                bottomVerticesBorder = sharedVertices(apicoBasalVertices(:, 3) < mean(midZ));
                
                apicalEdges = all(ismember(currentEdgesOfCell, upperVerticesBorder), 2);
                basalEdges = all(ismember(currentEdgesOfCell, bottomVerticesBorder), 2);
                lateralEdges = any(ismember(currentEdgesOfCell, upperVerticesBorder), 2) + any(ismember(currentEdgesOfCell, bottomVerticesBorder), 2) == 2;
                
                
                %% Get all apical and basal vertices
                upperZMinimum = (mean(midZ) + mean(Y.DataRow(upperVerticesBorder, 3))/10);
                bottomZMinimum = (mean(midZ) - mean(Y.DataRow(upperVerticesBorder, 3))/10);       

                obj.ApicalVertices{numCell} = vertcat(uniqueCurrentVertices(Y.DataRow(uniqueCurrentVertices, 3) > upperZMinimum), - find(obj.FaceCentres.DataRow(1:obj.FaceCentres.n, 3) > upperZMinimum));
                obj.ApicalBorderVertices{numCell} = ismember(obj.ApicalVertices{numCell}, upperVerticesBorder);
                
                obj.BasalVertices{numCell} = vertcat(uniqueCurrentVertices(Y.DataRow(uniqueCurrentVertices, 3) < bottomZMinimum), - find(obj.FaceCentres.DataRow(1:obj.FaceCentres.n, 3) < bottomZMinimum));
                obj.BasalBorderVertices{numCell} = ismember(obj.BasalVertices{numCell}, bottomVerticesBorder);
                
                obj.SubstrateForce{numCell} = zeros(length(obj.BasalVertices{numCell}), 1);
                
                %% Get apical and basal border vertices
                if sum(apicalEdges) > size(upperVerticesBorder, 1)
                    apicalVerticesIds = currentEdgesOfCell(apicalEdges, :);
                    [numElements, elements] = hist(apicalVerticesIds(:), unique(currentEdgesOfCell(apicalEdges, :)));
                    edgesToRemove = apicalVerticesIds(all(ismember(apicalVerticesIds, elements(numElements>2)), 2), :);
                    apicalEdges(ismember(currentEdgesOfCell, edgesToRemove, 'row')) = 0;
%                     apicalPixels = apicoBasalVertices(apicoBasalVertices(:, 3) > mean(midZ), :);
%                     %[newEdges] = boundaryOfCell(apicalPixels(:, 1:2));
%                     [newEdges] = boundaryOfCell(cmdscale(pdist(apicalPixels), 2));
%                     apicalEdges = ismember(sort(currentEdgesOfCell, 2), sort(apicalVertices(newEdges), 2), 'rows');
%                     
%                     if sum(apicalEdges) ~= size(apicalVertices, 1)
%                         error('CellClass:boundary issue');
%                     end
                end
                
                % 2:Basal 3:Apical 1:Lateral
                obj.EdgeLocation{numCell} = lateralEdges + 3*apicalEdges + 2*basalEdges;
                
            end
        end
    end
end