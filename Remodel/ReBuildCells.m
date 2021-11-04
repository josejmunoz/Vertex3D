function [Cell,nC, SCn, flag]=ReBuildCells(Cell,T,Y,X,SCn)
%% This function rebuilds Cells and Faces after doing local transformation (e.g 23flip.... )
% Loop i over cell
%   Loop j over the segment connecting the node (cell centre) i with neighbouring nodes (cells) j,  in the same time each segment (ij) correspond to a face shared by cell i and j.
%       Loop k over the vertices of face ij, which they correspond to the nodal tetrahedrons that go around the segment ij.
%        Then, face centres (Cell.FaceCentres) are placed in the centre of faces, except faces with three vertices which they do not have centres, since they are already triangles.


TetIDs=1:size(T.DataRow,1);
nC=[];
flag=false;
includedX=false(size(Cell.Int, 2),1);
for numCell = Cell.Int(ismember(Cell.Int,Cell.AssembleNodes))
    includedX(numCell)=true;
    Copy_cNodes=Cell.cNodes{numCell};
    Copy_Surface=Cell.Faces{numCell};
    Cell.nTotalTris=Cell.nTotalTris-size(Cell.Tris{numCell},1);
    
    %% Build Tet
    currentTets = any(ismember(T.DataRow,numCell),2);
    Cell.cTet{numCell}=T.DataRow(currentTets,:);
    Cell.cTetID{numCell} = TetIDs(currentTets);
    % Fill connected nodes
    Cell.cNodes{numCell} = unique(Cell.cTet{numCell});
    % Remove itself from the connections
    Cell.cNodes{numCell}(Cell.cNodes{numCell} == numCell) = [];
    
    %% Initiate cellular database
    nFaces=length(Cell.cNodes{numCell});
    Cell.Faces{numCell}.nFaces=nFaces;
    Cell.Faces{numCell}.FaceCentresID=zeros(nFaces,1);
    Cell.Faces{numCell}.Vertices=cell(nFaces,1);
    Cell.Faces{numCell}.Tris=cell(nFaces,1);
    Cell.Tris{numCell}=zeros(max(cellfun(@length, Cell.Tris)),3);
    Cell.Cv{numCell}=zeros(max(cellfun(@length, Cell.Cv)),2);
    nSurfTris=1;  % counter for cellular-number of SurfsTris
    nVertexBars=1;  % counter for cellular-number of vertex-bars
    
    %% Loop over all the faces of the cell
    for numNode = 1:length(Cell.cNodes{numCell})
        
        SurfAxes=[Cell.Int(numCell) Cell.cNodes{numCell}(numNode)];                               % line crossing the surface
        SurfTet=Cell.cTet{numCell}(sum(ismember(Cell.cTet{numCell},SurfAxes), 2) == 2, :);
        SurfTetID=Cell.cTetID{numCell}(sum(ismember(Cell.cTet{numCell},SurfAxes), 2) == 2);
        
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
        if length(remainingVertices)<1
            flag=true;
            return
        end
        SurfVertices(2)=remainingSurfTetID(remainingVertices(1));
        %remove the Found
        previousSurfTet=remainingSurfTet(remainingVertices(1),:);
        remainingSurfTet(remainingVertices(1),:)=[];
        remainingSurfTetID(remainingVertices(1))=[];
        for k=3:length(SurfVertices) %loop over Surf-vertices
            % Find the next, as the one who share three nodes with previous
            NextVertexlogic=sum(ismember(remainingSurfTet,previousSurfTet),2)==3;
            if length(remainingSurfTetID(NextVertexlogic))~=1
                flag=true;
                return
            end
            SurfVertices(k)=remainingSurfTetID(NextVertexlogic);
            %remove the Found
            previousSurfTet=remainingSurfTet(NextVertexlogic,:);
            remainingSurfTet(NextVertexlogic,:)=[];
            remainingSurfTetID(NextVertexlogic)=[];
        end
        
        %---build the center of the face (interpolation)
        % look for the center ID in the previous structure
        % (by cheking the corresponding Nodal connectivity )
        extrapolateFaceCentre = 0;
        numFaceCentresFaces = -1; % Change this to be the number of faceCentres faces (maybe Faces.n or idOppossedNode)

        % check if the centre i already built
        oppNode = Copy_cNodes==SurfAxes(2);
        
        idOppossedNode=Copy_Surface.FaceCentresID(oppNode);
        
        opposedVertices=[];
        if isempty(idOppossedNode) && ~isempty(nC)
            % if it is not found check if it is already added by another cell
            for numFace=1:length(nC)
                if all(ismember(SurfAxes,Cell.AllFaces.Nodes(nC(numFace),:)))
                    idOppossedNode=nC(numFace);
                    opposedVertices=Cell.AllFaces.Vertices{idOppossedNode};
                end
            end
        elseif ~isempty(idOppossedNode)
            opposedVertices=Copy_Surface.Vertices{Copy_cNodes==SurfAxes(2)};
        end
        
        if ~isempty(idOppossedNode) && length(opposedVertices)==1
            faceCentrePos=Cell.FaceCentres.DataRow(idOppossedNode,:);
        else
            faceCentrePos=sum(Y.DataRow(SurfVertices,:),1)/length(SurfVertices);
            if sum(ismember(SurfAxes,Cell.Int))==1 && extrapolateFaceCentre
                dir=(faceCentrePos-X(Cell.Int(numCell),:)); dir=dir/norm(dir);
                faceCentrePos=X(Cell.Int(numCell),:)+H.*dir;
            end
            %previousSurfTet=numFaceCentresFaces;
        end
        
        % % check orientation
        % v1=Y.DataRow(SurfVertices(1),:)-aux;
        % v2=Y.DataRow(SurfVertices(2),:)-aux;
        % if dot(cross(v1,v2),aux-X(Cell.Int(i),:))<0
        %     SurfVertices=flip(SurfVertices);
        % end
        Order=0;
        for iii=1:length(SurfVertices)
            if iii==length(SurfVertices)
                v1=Y.DataRow(SurfVertices(iii),:)-faceCentrePos;
                v2=Y.DataRow(SurfVertices(1),:)-faceCentrePos;
                Order=Order+dot(cross(v1,v2),faceCentrePos-X(Cell.Int(numCell),:))/length(SurfVertices);
            else
                v1=Y.DataRow(SurfVertices(iii),:)-faceCentrePos;
                v2=Y.DataRow(SurfVertices(iii+1),:)-faceCentrePos;
                Order=Order+dot(cross(v1,v2),faceCentrePos-X(Cell.Int(numCell),:))/length(SurfVertices);
            end
        end
        if Order<0
            SurfVertices=flip(SurfVertices);
        end
        
        % save face-data
        if ~isempty(idOppossedNode)
            Cell.AllFaces.Vertices{idOppossedNode}=SurfVertices;
            Cell.FaceCentres.DataRow(idOppossedNode,:)=faceCentrePos;
        else
            [Cell.FaceCentres,nCC]=Cell.FaceCentres.Add(faceCentrePos);
            nC=[nC nCC];
            SCn=SCn.Add(faceCentrePos);
            [Cell.AllFaces,idOppossedNode]=Cell.AllFaces.Add(SurfAxes,SurfVertices,Y.DataRow,Cell.FaceCentres.DataRow);
        end
        
        % Save surface vertices
        Cell.Faces{numCell}.Vertices{numNode}=SurfVertices;
        if length(SurfVertices)==3
            Cell.Tris{numCell}(nSurfTris,:)=[SurfVertices(1) SurfVertices(2) -SurfVertices(3)];
            Cell.Faces{numCell}.Tris{numNode}=[SurfVertices(1) SurfVertices(2) -SurfVertices(3)];
            auxCv=[SurfVertices(1) SurfVertices(2);
                SurfVertices(2) SurfVertices(3);
                SurfVertices(3) SurfVertices(1)];
            nSurfTris=nSurfTris+1;
            
            auxCv(ismember(auxCv,Cell.Cv{numCell},'rows') | ismember(flip(auxCv,2),Cell.Cv{numCell},'rows'),:)=[];
            Cell.Cv{numCell}(nVertexBars:nVertexBars+size(auxCv,1)-1,:)=auxCv;
            nVertexBars=nVertexBars+size(auxCv,1);
            
            % save Surface-center
            Cell.Faces{numCell}.FaceCentresID(numNode)=idOppossedNode;
        else
            % Build Tringles and CellCv
            auxCv=zeros(length(SurfVertices),2);
            Cell.Faces{numCell}.Tris{numNode}=zeros(length(SurfVertices),3);
            for h=2:length(SurfVertices)
                Cell.Tris{numCell}(nSurfTris,:)=[SurfVertices(h-1) SurfVertices(h) idOppossedNode];
                Cell.Faces{numCell}.Tris{numNode}(h-1,:)=[SurfVertices(h-1) SurfVertices(h) idOppossedNode];
                auxCv(2*(h-1)-1,:)=[SurfVertices(h-1) SurfVertices(h)];
                auxCv(2*(h-1),:)=[SurfVertices(h-1) -idOppossedNode];
                %             Cell.Cv{i}(count4)=
                nSurfTris=nSurfTris+1;
            end
            Cell.Tris{numCell}(nSurfTris,:)=[SurfVertices(end) SurfVertices(1) idOppossedNode];
            Cell.Faces{numCell}.Tris{numNode}(h,:)=[SurfVertices(end) SurfVertices(1) idOppossedNode];
            
            auxCv(2*h-1,:)=[SurfVertices(end) SurfVertices(1)];
            auxCv(2*h,:)=[SurfVertices(end) -idOppossedNode];
            % remove do duplicated bars
            auxCv(ismember(auxCv,Cell.Cv{numCell},'rows') | ismember(flip(auxCv,2),Cell.Cv{numCell},'rows'),:)=[];
            Cell.Cv{numCell}(nVertexBars:nVertexBars+size(auxCv,1)-1,:)=auxCv;
            %          count2=count2+1;
            nSurfTris=nSurfTris+1;
            nVertexBars=nVertexBars+size(auxCv,1);
            
            
            % save Surface-center
            Cell.Faces{numCell}.FaceCentresID(numNode)=idOppossedNode;
        end
    end
    Cell.Tris{numCell}(nSurfTris:end,:)=[];
    Cell.Cv{numCell}(nVertexBars:end,:)=[];
    Cell.nTotalTris=Cell.nTotalTris+size(Cell.Tris{numCell},1);
    
end
end
