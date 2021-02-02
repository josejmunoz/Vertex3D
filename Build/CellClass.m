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
        GhostCells            % Ghost cells
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
                Cell.GhostCells=false(nC, 1);
                Cell.ContractileForces=cell(nC, 1);
            end
        end
        
        %% Ablate cells
        % We assume there will only be one ablation. Thus, we could remove
        % the IDs from the middle.
        % TODO: Check if Cell.Int need to be in order in consecutive numbers
        function Cell = AblateCells(obj, cellsToRemove)
            obj.GhostCells(obj.Int == cellsToRemove) = true;
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
            obj.GhostCells(cellsToRemove) = [];
            
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
        
        function [obj, featuresTable] = exportTableWithCellFeatures(obj, Y)
            resolutionOfImage = 0.001;
            featuresTable = [];
            for numCell = obj.Int
                allVertices = [Y.DataRow; obj.FaceCentres.DataRow];
                
                verticesConnectionsOfCell = obj.Tris{numCell};
                verticesConnectionsOfCell(verticesConnectionsOfCell(:, 3)>=0, 3) = verticesConnectionsOfCell(verticesConnectionsOfCell(:, 3)>=0, 3) + size(allVertices, 1);
                verticesConnectionsOfCell(verticesConnectionsOfCell(:, 3)<0, 3) = abs(verticesConnectionsOfCell(verticesConnectionsOfCell(:, 3)<0, 3));
                %triangulationsOfCell = triangulation(verticesConnectionsOfCell, allVertices);
                
                DT = delaunayTriangulation(allVertices(unique(verticesConnectionsOfCell), :));
                [X,Y,Z] = meshgrid(min(allVertices(:)):resolutionOfImage:max(allVertices(:)));
                SI = pointLocation(DT,X(:),Y(:),Z(:));       %index of simplex (returns NaN for all points outside the convex hull)
                mask = ~isnan(SI); %binary
                mask = reshape(mask,size(X));
                currentTable = regionprops3(mask, 'all');
                currentTable.ID = numCell;
                featuresTable = [featuresTable; currentTable];
            end
        end
    end
end