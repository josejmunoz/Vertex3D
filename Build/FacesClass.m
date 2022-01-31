classdef FacesClass
    % Face Class->  this class stores the information of the faces of cells
    % and preforms face-based calculations  (e.g. Area, Energy). The main
    % objective of this class is to have the ability to remove and add
    % faces without the need to change the order and the index of all faces.
    % The energy --> WBexp = exp( Set.lambdaB*  ( 1 - Set.Beta*At/Set.BarrierTri0 )  );
    properties
        Energy         % - The energy of faces (array-structure ,  Size=[NumFaces 1]):
        %                   (e.g. Energy(i) -> The energy of the face i)
        %---------------------------------------------------------------------
        EnergyTri      % - The energy of face triangles (cell-structure ,  Size={NumFaces 1}):
        %                  (e.g. EnergyTri{i} -> array with the energy of the triangles of face i)
        %---------------------------------------------------------------------
        Vertices       % - The vertices of faces (cell-structure ,  Size={NumFaces 1}):
        %                  (e.g. Vertices{i} -> array with the vertices of face i
        %---------------------------------------------------------------------
        Area           % - The area of faces (array-structure ,  Size=[NumFaces 1]):
        %                  (e.g. Area(i) -> The area of face i)
        %---------------------------------------------------------------------
        AreaTri        % - The area of face triangles (cell-structure ,  Size={NumFaces 1}):
        %                   (e.g. AreaTri{i} -> array with the area of the triangles of face i)
        %---------------------------------------------------------------------
        Nodes          % - Face Nodes  (array-structure ,  Size=[NumFaces 2]):
        %                   (e.g. Nodes(i,:)=[n1 n2] -> are the two nodes which correspond to face i)
        %---------------------------------------------------------------------
        InterfaceType   % - Face Type  (array-structure ,  Size=[NumFaces 1]):
        %                    InterfaceType(i)=0 --> Face (i) is External
        %                    InterfaceType(i)=1 --> Face (i) is Internal (cell-cell interface)
        %                    InterfaceType(i)=2 --> Face (i) Substrate face (cell-substrate interface)
        %                    InterfaceType(i)=3 --> Face (i) is Cell-DebrisCell (cell-ablated cell interface)
        %---------------------------------------------------------------------
        NotEmpty        % - Logical array of size [NumFaces 1]
        %                    NotEmpty(i)= false -> Face (i) is empty,
        %                    NotEmpty(i)= true ->  Face (i) is not empty
        %---------------------------------------------------------------------
        nE              % - The Number of empty faces (storage spots)
        %---------------------------------------------------------------------
        EmptyList       % - List of empty Faces
        %---------------------------------------------------------------------
        n               % counter on the list of faces
    end
    
    %%
    methods
        function Array= FacesClass(S1)
            Array.Vertices=cell(S1,1);
            Array.AreaTri=cell(S1,1);
            Array.Area=zeros(S1,1);
            Array.Nodes=zeros(S1,2);
            Array.Energy=zeros(S1,1);
            Array.EnergyTri=cell(S1,1);
            Array.NotEmpty=false(S1,1);
            Array.InterfaceType=zeros(S1,1);
            Array.EmptyList=nan(S1,1);
            Array.n=0;
            Array.nE=0;
        end
        
        %------------- Get the actual number of faces ---------------------
        function s=Size(obj)
            s=obj.n-obj.nE;
        end
        
        %
        %------------- Add new face ---------------------------------------
        function [obj,list] = Add(obj,Nodes,Vertices,Y,SurfsCenters)
            if obj.nE>0
                % There are empty faces (storage spots) -> refill them first
                indexToAdd = obj.EmptyList(obj.nE);
            else
                % There are no empty faces (storage spots) -> Add the new face at the end of the list
                indexToAdd = obj.n+1;
            end
            
            obj.Vertices{indexToAdd}=Vertices;
            obj.Nodes(indexToAdd,:)=Nodes;
            obj.Energy(indexToAdd,:)=0;
            
            % Compute Area
            obj.AreaTri{indexToAdd}=zeros(length(Vertices),1);
            obj.EnergyTri{indexToAdd}=zeros(length(Vertices),1);
            for j=1:length(Vertices)-1
                YTri=[Y(Vertices([j j+1]),:); SurfsCenters(indexToAdd,:)];
                obj.AreaTri{indexToAdd}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                %                         obj.AreaTri0{indexToAdd}(j)=obj.AreaTri{indexToAdd}(j);
            end
            YTri=[Y(Vertices([j+1 1]),:); SurfsCenters(indexToAdd,:)];
            obj.AreaTri{indexToAdd}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
            obj.Area(indexToAdd)=sum(obj.AreaTri{indexToAdd});
            
            if obj.nE>0
                % There are empty faces (storage spots) -> refill them first
                list=obj.EmptyList(obj.nE);
                obj.NotEmpty(obj.EmptyList(obj.nE))=true;
                obj.EmptyList(obj.nE)=nan;
                obj.nE=obj.nE-1;
            else
                % There are no empty faces (storage spots) -> Add the new face at the end of the list
                list=obj.n+1;
                obj.NotEmpty(obj.n+1)=true;
                obj.n=obj.n+1;
            end
        end
        
        %-------------Remove face -------------------------------------------
        function [obj] = Remove(obj,V)
            for i=1:length(V)
                obj.Vertices{V(i)}=[];
            end
            obj.Nodes(V,:)=nan;
            obj.Energy(V)=nan;
            obj.Area(V)=nan;
            obj.InterfaceType(V)=0;
            obj.AreaTri{V}=[];
            obj.NotEmpty(V)=false;
            obj.EmptyList(obj.nE+1:obj.nE+length(V))=V;
            obj.nE=obj.nE+length(V);
        end
        
        %-------------Compute the area of faces ----------------------------
        function [obj]=ComputeAreaTri(obj,Y,SurfsCenters)
            for i=1:obj.n
                if obj.NotEmpty(i)
                    obj.AreaTri{i}=zeros(length(obj.Vertices{i}),1);
                    for j=1:length(obj.Vertices{i})-1
                        YTri=[Y(obj.Vertices{i}([j j+1]),:); SurfsCenters(i,:)];
                        obj.AreaTri{i}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    end
                    YTri=[Y(obj.Vertices{i}([j+1 1]),:); SurfsCenters(i,:)];
                    obj.AreaTri{i}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    obj.Area(i)=sum(obj.AreaTri{i});
                else
                    obj.AreaTri{i}=[];
                end
            end
        end
        
        
        %-------------Compute the energy barrier of faces -----------------
        % The energy --> WBexp = exp( Set.lambdaB*  ( 1 - Set.Beta*At/Set.BarrierTri0 )  );
        function [obj]=ComputeEnergy(obj,Set)
            for i=1:obj.n
                if obj.NotEmpty(i)
                    obj.EnergyTri{i}=zeros(size(obj.Vertices{i}));
                    for j=1:length(obj.Vertices{i})
                        obj.EnergyTri{i}(j)= exp(Set.lambdaB*(1-Set.Beta*obj.AreaTri{i}(j)/Set.BarrierTri0));
                    end
                    obj.Energy(i)=sum(obj.EnergyTri{i});
                else
                    obj.EnergyTri{i}=[];
                    obj.Energy(i)=0;
                end
            end
        end
        
        %-------------Obtain the type of Faces -----------------
        function [obj]=CheckInteriorFaces(obj, Cell)
            for numFace=1:obj.n
                if obj.NotEmpty(numFace)
                    if all(ismember(obj.Nodes(numFace,:),Cell.Int))
                        % cell-cell face
                        obj.InterfaceType(numFace)=1;
                    else
                        if Cell.FaceCentres(numFace, 3) >= 0
                            % external/apical face
                            obj.InterfaceType(numFace)=0;
                        else
                            % basal face
                            obj.InterfaceType(numFace)=2;
                        end
                    end
                end
            end
        end
    end
end