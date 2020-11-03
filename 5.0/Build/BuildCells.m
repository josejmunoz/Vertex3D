function [Cv,Cell,SharedFaces]=BuildCells(TT,Y,X,xInternal,xExternal,H)
nC=length(xInternal);
Cell = CellClass(X,nC,xInternal,xExternal);

% -------
TetIDs=1:size(TT,1);
%% Initiate global database
Cv=zeros(Cell.n*16*8,2);
countVertexBarElements=1; % counter for the number of vertex-bar Elemets 

Cell.SurfsCenters=zeros(size(X,1)*4,3); % the posiion of Cell-surface centers
SharedFaces = FacesClass(size(X,1)*4);
countRatioSurfsCenters_Faces=1; % counter for the number of SurfsCenters / faces; 
countTotalSurfsTris=0; % counter for total number of SurfsTris 

Includedx=false(size(Cell.Int,2),1);
%% Build Interior Cells
for interiorCell=1:length(Cell.Int) % loop over  Interior cells
    Includedx(interiorCell)=true;
    % i Should be Cell.Int(i) when boundary nodes are included
    
    % ----- Build Tet
%     Cell.cTris{i}=[];
    Cell.cTet{interiorCell}=TT(any(ismember(TT,Cell.Int(interiorCell)),2),:);
    Cell.cTetID{interiorCell}=TetIDs(any(ismember(TT,Cell.Int(interiorCell)),2));
    Cell.cNodes{interiorCell}=unique(Cell.cTet{interiorCell}); 
    Cell.cNodes{interiorCell}(Cell.cNodes{interiorCell}==Cell.Int(interiorCell))=[];
    
    % ----- Build Surfaces and Trinagles   
    %-- Initiate cellular database
%    Cell.CellSurfsCentersID{i}=zeros(length(Cell.cNodes{i}),1);
    nS=length(Cell.cNodes{interiorCell});
    Cell.Surfaces{interiorCell}.nSurfaces=nS;
    Cell.Surfaces{interiorCell}.SurfaceCentersID=zeros(nS,1);
    Cell.Surfaces{interiorCell}.SurfaceVertices=cell(nS,1);
    Cell.Surfaces{interiorCell}.Tris=cell(nS,1);

    Cell.Tris{interiorCell}=zeros(120,3);
    countCell_SurfsTris=1;  % counter for cellular-number of SurfsTris   
    Cell.Cv{interiorCell}=zeros(16*8,2);
    countCell_VertexBars=1;  % counter for cellular-number of vertex-bars  
    
    for numCellSurface=1:length(Cell.cNodes{interiorCell}) % loop over cell-Surfaces
        
        SurfAxes=[Cell.Int(interiorCell) Cell.cNodes{interiorCell}(numCellSurface)];                               % line crossing the surface
        SurfTet=Cell.cTet{interiorCell}(sum(ismember(Cell.cTet{interiorCell},SurfAxes),2)==2,:);      
        SurfTetID=Cell.cTetID{interiorCell}(sum(ismember(Cell.cTet{interiorCell},SurfAxes),2)==2);
        
        %--- Order Surface Vertices\Tets as loop 
        SurfVertices=zeros(length(SurfTetID),1);   
        AuxSurfTet=SurfTet; 
        AuxSurfTetID=SurfTetID; % To be arranged
        % Take the first 
        SurfVertices(1)=SurfTetID(1);
        AuxSurfTet(1,:)=[];
        AuxSurfTetID(1)=[];
        % Take the nearby regardless of the direction
        NextVertexlogic=sum(ismember(AuxSurfTet,SurfTet(1,:)),2)==3;
        aux1=1:length(NextVertexlogic); aux1=aux1(NextVertexlogic);
        SurfVertices(2)=AuxSurfTetID(aux1(1));
        %remove the Found
        aux2=AuxSurfTet(aux1(1),:);
        AuxSurfTet(aux1(1),:)=[];
        AuxSurfTetID(aux1(1))=[];
        for numSurfVertex=3:length(SurfVertices) %loop over Surf-vertices
            % Find the next, as the one who share three nodes with previous
            NextVertexlogic=sum(ismember(AuxSurfTet,aux2),2)==3;
            SurfVertices(numSurfVertex)=AuxSurfTetID(NextVertexlogic);
            %remove the Found
            aux2=AuxSurfTet(NextVertexlogic,:);
            AuxSurfTet(NextVertexlogic,:)=[];
            AuxSurfTetID(NextVertexlogic)=[];
        end 
        
        % build the center of the face (interpolation)
          
         % check if the center i already built
         oppNode=Cell.Int==SurfAxes(2); 
         if Includedx(oppNode)
             cID=Cell.Surfaces{oppNode}.SurfaceCentersID(Cell.cNodes{oppNode}==SurfAxes(1));
             aux2=cID;
             aux=Cell.SurfsCenters(aux2,:);
         else 
             aux=sum(Y.DataRow(SurfVertices,:),1)/length(SurfVertices);
             if sum(ismember(SurfAxes,Cell.Int))==1
                 dir=(aux-X(Cell.Int(interiorCell),:)); dir=dir/norm(dir);
                 aux=X(Cell.Int(interiorCell),:)+H.*dir;
             end
              aux2=countRatioSurfsCenters_Faces;
         end 
         
         
         
        % check  Orientation
%         v1=Y.DataRow(SurfVertices(1),:)-aux;
%         v2=Y.DataRow(SurfVertices(2),:)-aux;
        if interiorCell==5
            malik=1;
        end 
            
        Order=0;
        for iii=1:length(SurfVertices)
            if iii==length(SurfVertices)
                v1=Y.DataRow(SurfVertices(iii),:)-aux;
                v2=Y.DataRow(SurfVertices(1),:)-aux;
                Order=Order+dot(cross(v1,v2),aux-X(Cell.Int(interiorCell),:))/length(SurfVertices);
            else 
                v1=Y.DataRow(SurfVertices(iii),:)-aux;
                v2=Y.DataRow(SurfVertices(iii+1),:)-aux;
                Order=Order+dot(cross(v1,v2),aux-X(Cell.Int(interiorCell),:))/length(SurfVertices);
            end 
        end 
        if Order<0
           SurfVertices=flip(SurfVertices);
        end 
        
        % Save surface vertices 
        Cell.Surfaces{interiorCell}.SurfaceVertices{numCellSurface}=SurfVertices;
        
         % Build Tringles and CellCv
         % Build Tringles and CellCv
         if length(SurfVertices)==3 % && false
             Cell.Tris{interiorCell}(countCell_SurfsTris,:)=[SurfVertices(1) SurfVertices(2) -SurfVertices(3)];
             Cell.Surfaces{interiorCell}.Tris{numCellSurface}=[SurfVertices(1) SurfVertices(2) -SurfVertices(3)];
             auxCv=[SurfVertices(1) SurfVertices(2);
                           SurfVertices(2) SurfVertices(3);
                           SurfVertices(3) SurfVertices(1)];
             countTotalSurfsTris=countTotalSurfsTris+1;
             countCell_SurfsTris=countCell_SurfsTris+1;
             
             auxCv(ismember(auxCv,Cell.Cv{interiorCell},'rows') | ismember(flip(auxCv,2),Cell.Cv{interiorCell},'rows'),:)=[];
             Cell.Cv{interiorCell}(countCell_VertexBars:countCell_VertexBars+size(auxCv,1)-1,:)=auxCv;
             countCell_VertexBars=countCell_VertexBars+size(auxCv,1);
             
            % save Surface-center
             if Includedx(oppNode)
                 Cell.Surfaces{interiorCell}.SurfaceCentersID(numCellSurface)=cID;
             else 
                 Cell.SurfsCenters(countRatioSurfsCenters_Faces,:)=aux;
                 Cell.Surfaces{interiorCell}.SurfaceCentersID(numCellSurface)=countRatioSurfsCenters_Faces;
                 % Save Face-data base 
                  SharedFaces = SharedFaces.Add(SurfAxes,SurfVertices,Y.DataOrdered,Cell.SurfsCenters);
                 countRatioSurfsCenters_Faces=countRatioSurfsCenters_Faces+1;
             end 
         else 
             auxCv=zeros(length(SurfVertices),2);
             Cell.Surfaces{interiorCell}.Tris{numCellSurface}=zeros(length(SurfVertices),3);
             for h=2:length(SurfVertices)
                Cell.Tris{interiorCell}(countCell_SurfsTris,:)=[SurfVertices(h-1) SurfVertices(h) aux2];
                Cell.Surfaces{interiorCell}.Tris{numCellSurface}(h-1,:)=[SurfVertices(h-1) SurfVertices(h) aux2];
                auxCv(h-1,:)=[SurfVertices(h-1) SurfVertices(h)];
                countTotalSurfsTris=countTotalSurfsTris+1;
                countCell_SurfsTris=countCell_SurfsTris+1;
             end 
             Cell.Tris{interiorCell}(countCell_SurfsTris,:)=[SurfVertices(end) SurfVertices(1) aux2];
             Cell.Surfaces{interiorCell}.Tris{numCellSurface}(h,:)=[SurfVertices(end) SurfVertices(1) aux2];

             auxCv(h,:)=[SurfVertices(end) SurfVertices(1)];
             % remove do duplicated bars
             auxCv(ismember(auxCv,Cell.Cv{interiorCell},'rows') | ismember(flip(auxCv,2),Cell.Cv{interiorCell},'rows'),:)=[];
             Cell.Cv{interiorCell}(countCell_VertexBars:countCell_VertexBars+size(auxCv,1)-1,:)=auxCv;
             countTotalSurfsTris=countTotalSurfsTris+1;
             countCell_SurfsTris=countCell_SurfsTris+1;
             countCell_VertexBars=countCell_VertexBars+size(auxCv,1);
              % save Surface-center
             if Includedx(oppNode)
                 Cell.Surfaces{interiorCell}.SurfaceCentersID(numCellSurface)=cID;
             else 
                 Cell.SurfsCenters(countRatioSurfsCenters_Faces,:)=aux;
                 Cell.Surfaces{interiorCell}.SurfaceCentersID(numCellSurface)=countRatioSurfsCenters_Faces;
                 % Save Face-data base 
                  SharedFaces = SharedFaces.Add(SurfAxes,SurfVertices,Y.DataOrdered,Cell.SurfsCenters);
                 countRatioSurfsCenters_Faces=countRatioSurfsCenters_Faces+1;
             end 
         end 

         

    end 
    Cell.Tris{interiorCell}(countCell_SurfsTris:end,:)=[];
    Cell.Cv{interiorCell}(countCell_VertexBars:end,:)=[];
    
    Cell.CvID{interiorCell}=countVertexBarElements:countVertexBarElements+size(Cell.Cv{interiorCell},1)-1;
    Cv(countVertexBarElements:countVertexBarElements+size(Cell.Cv{interiorCell},1)-1,:)=Cell.Cv{interiorCell};
    countVertexBarElements=countVertexBarElements+size(Cell.Cv{interiorCell},1);
end 
Cell.SurfsCenters(countRatioSurfsCenters_Faces:end,:)=[];
Cv(countVertexBarElements:end,:)=[];
Cell.nTotalTris=countTotalSurfsTris;



%% change type of data strucutre (should be done in the beginning)

aux=Cell.SurfsCenters;
Cell.SurfsCenters=DynamicArray(size(X,1)*8,3);
Cell.SurfsCenters=Cell.SurfsCenters.Add(aux);
[Cell]=BuildEdges(Cell,Y);
%% Compute Cells volume 
[Cell]=ComputeCellVolume(Cell,Y);
Cell.Vol0=Cell.Vol;
Cell.SArea0=Cell.SArea;
for interiorCell=1:Cell.n
%     Cell.SAreaTri0{i}=Cell.SAreaTri{i}*1e-2;
    Cell.SAreaTri0{interiorCell}=ones(size(Cell.SAreaTri{interiorCell}))*1e-3;
%     Cell.SAreaTri0{i}=Cell.SAreaTri{i};

end 
Cell.SAreaFace0=Cell.SAreaFace;

SharedFaces=SharedFaces.ComputeAreaTri(Y.DataRow,Cell.SurfsCenters.DataRow);


end 



%%



