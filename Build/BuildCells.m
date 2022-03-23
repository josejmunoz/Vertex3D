function [Cv,Cell]=BuildCells(T,Y,X,xInternal,H, extrapolateFaceCentre)
% Basically in this function  the cells are built by doing the following loops
% Loop i over cell
%   Loop j over the segment connecting the node (cell centre) i with neighbouring nodes (cells) j,  in the same time each segment (ij) correspond to a face shared by cell i and j.
%       Loop k over the vertices of face ij, which they correspond to the nodal tetrahedrons that go around the segment ij.
% Then, face centres (Cell.FaceCentres) are obviously placed in the centre of faces, except faces with three vertices which they do not have centres, since they are already triangles.



nC=length(xInternal);
Cell = CellClass(nC,xInternal);
% -------
TetIDs=1:size(T,1);
%% Initiate global database
Cv=zeros(Cell.n*16*8,2);
numVertexBarElem=1; % counter for the number of vertex-bar Elements
Cell.FaceCentres=zeros(size(X,1)*4,3); % the position of Cell-surface centres
Cell.AllFaces = FacesClass(size(X,1)*4);
numFaceCentresFaces=1; % counter for the number of FaceCentres / faces;
numTotalSurfsTris=0; % counter for total number of SurfsTris

Includedx=false(size(Cell.Int, 2),1);
%% Build Interior Cells
for i=1:length(Cell.Int) % loop over  Interior cells
    Includedx(i)=true;
    Cell.Centre(i, :) = X(i, :);
    % i Should be Cell.Int(i) when boundary nodes are included
    
    %% Build Tetrahedra of the nodes
    currentTetrahedra = any(ismember(T,Cell.Int(i)),2);
    Cell.cTet{i}=T(currentTetrahedra,:);
    Cell.cTetID{i}=TetIDs(currentTetrahedra);
    Cell.cNodes{i}=unique(Cell.cTet{i});
    Cell.cNodes{i}(Cell.cNodes{i}==Cell.Int(i))=[];
    
    %% Build Faces and Triangles
    nS=length(Cell.cNodes{i});
    Cell.Faces{i}.nFaces=nS;
    Cell.Faces{i}.FaceCentresID=zeros(nS,1);
    Cell.Faces{i}.Vertices=cell(nS,1);
    Cell.Faces{i}.Tris=cell(nS,1);
    
    Cell.Tris{i}=zeros(120,3);
    numCell_SurfsTris=1;  % counter for cellular-number of SurfsTris
    Cell.Cv{i}=zeros(16*8,2);
    numCell_VertexBars=1;  % counter for cellular-number of vertex-bars
    
    %% Cell-faces
    for j=1:length(Cell.cNodes{i})
        % line crossing the surface
        SurfAxes=[Cell.Int(i) Cell.cNodes{i}(j)];
        SurfTet=Cell.cTet{i}(sum(ismember(Cell.cTet{i},SurfAxes),2)==2,:);
        SurfTetID=Cell.cTetID{i}(sum(ismember(Cell.cTet{i},SurfAxes),2)==2);
        
        %% Order Surface Vertices\Tets as loop
        SurfVertices=zeros(length(SurfTetID),1);
        remainingSurfTet=SurfTet;
        remainingSurfTetID=SurfTetID; % To be arranged
        % Take the first
        SurfVertices(1)=SurfTetID(1);
        remainingSurfTet(1,:)=[];
        remainingSurfTetID(1)=[];
        % Take the nearby regardless of the direction
        NextVertexlogic=sum(ismember(remainingSurfTet,SurfTet(1,:)),2)==3;
        remainingVertices=1:length(NextVertexlogic); 
        remainingVertices=remainingVertices(NextVertexlogic);
        SurfVertices(2)=remainingSurfTetID(remainingVertices(1));
        %remove the Found
        previousSurfTet=remainingSurfTet(remainingVertices(1),:);
        remainingSurfTet(remainingVertices(1),:)=[];
        remainingSurfTetID(remainingVertices(1))=[];
        for k=3:length(SurfVertices) %loop over Surf-vertices
            % Find the next, as the one who share three nodes with previous
            NextVertexlogic=sum(ismember(remainingSurfTet,previousSurfTet),2)==3;
            SurfVertices(k)=remainingSurfTetID(NextVertexlogic);
            %remove the Found
            previousSurfTet=remainingSurfTet(NextVertexlogic,:);
            remainingSurfTet(NextVertexlogic,:)=[];
            remainingSurfTetID(NextVertexlogic)=[];
        end
        
        %% Build the centre of the face (interpolation)
        [SurfVertices, previousSurfTet, faceCentrePos, oppNode, cID] = BuildFaceCentre(i, SurfAxes, Cell, X, Y, SurfVertices, Includedx, H, numFaceCentresFaces, extrapolateFaceCentre);
        
        % Save Face vertices
        Cell.Faces{i}.Vertices{j}=SurfVertices;
        
        % save Face-centre
        if Includedx(oppNode)
            Cell.Faces{i}.FaceCentresID(j)=cID;
        else
            Cell.FaceCentres(numFaceCentresFaces,:)=faceCentrePos;
            Cell.Faces{i}.FaceCentresID(j)=numFaceCentresFaces;
            % Save Face-data base
            Cell.AllFaces = Cell.AllFaces.Add(SurfAxes,SurfVertices,Y.DataOrdered,Cell.FaceCentres);
            numFaceCentresFaces=numFaceCentresFaces+1;
        end
        
        
        % Build Triangles and CellCv
        if length(SurfVertices)==3 % && false
            % Build Triangles and CellCv
            Cell.Tris{i}(numCell_SurfsTris,:)=[SurfVertices(1) SurfVertices(2) -SurfVertices(3)];
            Cell.Faces{i}.Tris{j}=[SurfVertices(1) SurfVertices(2) -SurfVertices(3)];
            auxCv=[SurfVertices(1) SurfVertices(2);
                SurfVertices(2) SurfVertices(3);
                SurfVertices(3) SurfVertices(1)];
            numTotalSurfsTris=numTotalSurfsTris+1;
            numCell_SurfsTris=numCell_SurfsTris+1;
            
            auxCv(ismember(auxCv,Cell.Cv{i},'rows') | ismember(flip(auxCv,2),Cell.Cv{i},'rows'),:)=[];
            Cell.Cv{i}(numCell_VertexBars:numCell_VertexBars+size(auxCv,1)-1,:)=auxCv;
            numCell_VertexBars=numCell_VertexBars+size(auxCv,1);
        else
            auxCv=zeros(length(SurfVertices)*2,2);
            Cell.Faces{i}.Tris{j}=zeros(length(SurfVertices),3);
            for h=2:length(SurfVertices)
                Cell.Tris{i}(numCell_SurfsTris,:)=[SurfVertices(h-1) SurfVertices(h) previousSurfTet];
                Cell.Faces{i}.Tris{j}(h-1,:)=[SurfVertices(h-1) SurfVertices(h) previousSurfTet];
                auxCv(2*(h-1)-1,:)=[SurfVertices(h-1) SurfVertices(h)];
                auxCv(2*(h-1),:)=[SurfVertices(h-1) -previousSurfTet];
                numTotalSurfsTris=numTotalSurfsTris+1;
                numCell_SurfsTris=numCell_SurfsTris+1;
            end
            Cell.Tris{i}(numCell_SurfsTris,:)=[SurfVertices(end) SurfVertices(1) previousSurfTet];
            Cell.Faces{i}.Tris{j}(h,:)=[SurfVertices(end) SurfVertices(1) previousSurfTet];
            
            auxCv(2*h-1,:)=[SurfVertices(end) SurfVertices(1)];
            auxCv(2*h,:)=[SurfVertices(end) -previousSurfTet];
            % remove do duplicated bars
            auxCv(ismember(auxCv,Cell.Cv{i},'rows') | ismember(flip(auxCv,2),Cell.Cv{i},'rows'),:)=[];
            Cell.Cv{i}(numCell_VertexBars:numCell_VertexBars+size(auxCv,1)-1,:)=auxCv;
            numTotalSurfsTris=numTotalSurfsTris+1;
            numCell_SurfsTris=numCell_SurfsTris+1;
            numCell_VertexBars=numCell_VertexBars+size(auxCv,1);
        end
    end
    Cell.Tris{i}(numCell_SurfsTris:end,:)=[];
    Cell.Cv{i}(numCell_VertexBars:end,:)=[];
    
    Cell.CvID{i}=numVertexBarElem:numVertexBarElem+size(Cell.Cv{i},1)-1;
    Cv(numVertexBarElem:numVertexBarElem+size(Cell.Cv{i},1)-1,:)=Cell.Cv{i};
    numVertexBarElem=numVertexBarElem+size(Cell.Cv{i},1);
end
Cell.FaceCentres(numFaceCentresFaces:end,:)=[];
Cv(numVertexBarElem:end,:)=[];
Cell.nTotalTris=numTotalSurfsTris;
Cell.Centre0 = Cell.Centre;

%% change type of data structure (should be done in the beginning)

faceCentrePos=Cell.FaceCentres;
Cell.FaceCentres=DynamicArray(size(X,1)*8,3);
Cell.FaceCentres=Cell.FaceCentres.Add(faceCentrePos);
[Cell]=BuildEdges(Cell,Y);

Cell.FaceCentres0 = Cell.FaceCentres;
Cell.Tris0 = Cell.Tris;

%% Compute Cells volume
[Cell]=ComputeCellVolume(Cell,Y);
Cell.Vol0=Cell.Vol;
% volumeDevMean = Cell.Vol - mean(Cell.Vol);
% Cell.Vol0 = Cell.Vol - (volumeDevMean/1.5);
Cell.SArea0=Cell.SArea;
for i=1:Cell.n
    Cell.SAreaTrin{i}=Cell.SAreaTri{i};
end
Cell.SAreaFace0=Cell.SAreaFace;
Cell.AllFaces=Cell.AllFaces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.ComputePerimeterTri(Y.DataRow,Cell.FaceCentres.DataRow);

%% Apico-basal distinction
[Cell] = Cell.computeEdgeLocation(Y);

%% Contractility L_0: Only when builiding the cells to avoid issues in remodelling
[Cell, uniqueEdges] = Cell.computeEdgeLengths(Y);
allEdges = vertcat(Cell.EdgeLengths{:});
allLocations = vertcat(Cell.EdgeLocation{:});
Cell.EdgeLengths0_average = mean(allEdges(allLocations == 3 | allLocations == 2));
Cell.EdgeLengths0_lateralAverage = mean(allEdges(allLocations == 1))/2;
Cell.EdgeLengthsn = Cell.EdgeLengths;

end