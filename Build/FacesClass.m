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
        %---------------------------------------------------------------------
        NotEmpty        % - Logical array of size [NumFaces 1]
        %                    NotEmpty(i)= false -> Face (i) is empty,
        %                    NotEmpty(i)= true ->  Face (i) is not empty
        %---------------------------------------------------------------------
        nE              % - The Number of empty faces (storage spots)
        %---------------------------------------------------------------------
        EmptyList       % - List of empty Faces
        %---------------------------------------------------------------------
        V3              % - Faces with 3 vertices ( Logical array of size [NumFaces 1] )
        %                     V3(i)= true -> Face (i) has  3 vertices.
        %                     V3(i)= false -> Face (i) dose not has 3 vertices.
        %---------------------------------------------------------------------
        V4              % - Faces with 4 vertices ( Logical array of size [NumFaces 1] )
        %                      V4(i)= true -> Face (i) has 4 vertices.
        %                      V4(i)= false -> Face (i) dose not has 4 vertices.
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
            Array.V3=false(S1,1);
            Array.V4=false(S1,1);
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
                obj.Vertices{obj.EmptyList(obj.nE)}=Vertices;
                obj.Nodes(obj.EmptyList(obj.nE),:)=Nodes;
                obj.Energy(obj.EmptyList(obj.nE),:)=0;
                if length(Vertices)==3
                    obj.V3(obj.EmptyList(obj.nE))=true;
                    obj.V4(obj.EmptyList(obj.nE))=false;
                    % Compute Area
                    YTri=Y(Vertices,:);
                    obj.AreaTri{obj.EmptyList(obj.nE)}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    obj.EnergyTri{obj.EmptyList(obj.nE)}=0;
                    obj.Area(obj.EmptyList(obj.nE))=sum(obj.AreaTri{obj.EmptyList(obj.nE)});
                    list=obj.EmptyList(obj.nE);
                    
                    obj.NotEmpty(obj.EmptyList(obj.nE))=true;
                    obj.EmptyList(obj.nE)=nan;
                    obj.nE=obj.nE-1;
                else
                    if length(Vertices)==4
                        obj.V4(obj.EmptyList(obj.nE))=true;
                    else
                        obj.V4(obj.EmptyList(obj.nE))=false;
                    end
                    obj.V3(obj.EmptyList(obj.nE))=false;
                    % Compute Area
                    obj.AreaTri{obj.EmptyList(obj.nE)}=zeros(length(Vertices),1);
                    obj.EnergyTri{obj.EmptyList(obj.nE)}=zeros(length(Vertices),1);
                    for j=1:length(Vertices)-1
                        YTri=[Y(Vertices([j j+1]),:); SurfsCenters(obj.EmptyList(obj.nE),:)];
                        obj.AreaTri{obj.EmptyList(obj.nE)}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                        %                         obj.AreaTri0{obj.EmptyList(obj.nE)}(j)=obj.AreaTri{obj.EmptyList(obj.nE)}(j);
                    end
                    YTri=[Y(Vertices([j+1 1]),:); SurfsCenters(obj.EmptyList(obj.nE),:)];
                    obj.AreaTri{obj.EmptyList(obj.nE)}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    obj.Area(obj.EmptyList(obj.nE))=sum(obj.AreaTri{obj.EmptyList(obj.nE)});
                                        
                    list=obj.EmptyList(obj.nE);
                    obj.NotEmpty(obj.EmptyList(obj.nE))=true;
                    obj.EmptyList(obj.nE)=nan;
                    obj.nE=obj.nE-1;
                end
                
            else
                % There are no empty faces (storage spots) -> Add the new face at the end of the list
                obj.Vertices{obj.n+1}=Vertices;
                obj.Nodes(obj.n+1,:)=Nodes;
                obj.Energy(obj.n+1)=0;
                if length(Vertices)==3
                    obj.V3(obj.n+1)=true;
                    obj.V4(obj.n+1)=false;
                    % Compute Area
                    YTri=Y(Vertices,:);
                    obj.AreaTri{obj.n+1}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    obj.EnergyTri{obj.n+1}=0;
                    obj.Area(obj.n+1)=sum(obj.AreaTri{obj.n+1});
                    list=obj.n+1;
                    obj.NotEmpty(obj.n+1)=true;
                    obj.n=obj.n+1;
                    
                else
                    obj.V3(obj.n+1)=false;
                    if length(Vertices)==4
                        obj.V4(obj.n+1)=true;
                    else
                        obj.V4(obj.n+1)=false;
                    end
                    % Compute Area
                    obj.AreaTri{obj.n+1}=zeros(length(Vertices),1);
                    obj.EnergyTri{obj.n+1}=zeros(length(Vertices),1);
                    for j=1:length(Vertices)-1
                        YTri=[Y(Vertices([j j+1]),:); SurfsCenters(obj.n+1,:)];
                        obj.AreaTri{obj.n+1}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    end
                    YTri=[Y(Vertices([j+1 1]),:); SurfsCenters(obj.n+1,:)];
                    obj.AreaTri{obj.n+1}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                    obj.Area(obj.n+1)=sum(obj.AreaTri{obj.n+1});
                    
                    list=obj.n+1;
                    obj.NotEmpty(obj.n+1)=true;
                    obj.n=obj.n+1;
                end
                
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
            obj.V3(V)=false;
            obj.V4(V)=false;
            obj.EmptyList(obj.nE+1:obj.nE+length(V))=V;
            obj.nE=obj.nE+length(V);
        end
        
        %-------------Compute the area of faces ----------------------------
        function [obj]=ComputeAreaTri(obj,Y,SurfsCenters)
            for i=1:obj.n
                if obj.NotEmpty(i)
                    if length(obj.Vertices{i})==3
                        YTri=Y(obj.Vertices{i},:);
                        obj.AreaTri{i}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                        obj.Area(i)=sum(obj.AreaTri{i});
                    else
                        obj.AreaTri{i}=zeros(length(obj.Vertices{i}),1);
                        for j=1:length(obj.Vertices{i})-1
                            YTri=[Y(obj.Vertices{i}([j j+1]),:); SurfsCenters(i,:)];
                            obj.AreaTri{i}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                        end
                        YTri=[Y(obj.Vertices{i}([j+1 1]),:); SurfsCenters(i,:)];
                        obj.AreaTri{i}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                        obj.Area(i)=sum(obj.AreaTri{i});
                    end
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
                    if length(obj.Vertices{i})==3
                        obj.EnergyTri{i}=exp(Set.lambdaB*(1-Set.Beta*obj.AreaTri{i}/Set.BarrierTri0));
                        obj.Energy(i)=sum(obj.EnergyTri{i});
                    else
                        obj.EnergyTri{i}=zeros(size(obj.Vertices{i}));
                        for j=1:length(obj.Vertices{i})
                            obj.EnergyTri{i}(j)= exp(Set.lambdaB*(1-Set.Beta*obj.AreaTri{i}(j)/Set.BarrierTri0));
                        end
                        obj.Energy(i)=sum(obj.EnergyTri{i});
                    end
                else
                    obj.EnergyTri{i}=[];
                    obj.Energy(i)=0;
                end
            end
        end
        
        %-------------Obtain the type of Faces -----------------
        function [obj]=CheckInteriorFaces(obj,XgID,XgSub)
            if nargin==2
                for i=1:obj.n
                    if obj.NotEmpty(i)
                        if any(ismember(obj.Nodes(i,:),XgID))
                            % external face
                            obj.InterfaceType(i)=0;
                        else
                            % cell-cell face
                            obj.InterfaceType(i)=1;
                        end
                    end
                end
            else
                XgID(XgID==XgSub)=[];
                for i=1:obj.n
                    if obj.NotEmpty(i)
                        if any(ismember(obj.Nodes(i,:),XgID))
                            % external face
                            obj.InterfaceType(i)=0;
                        elseif any(ismember(obj.Nodes(i,:),XgSub))
                            % cell-substrate face
                            obj.InterfaceType(i)=2;
                        else
                            % cell-cell face
                            obj.InterfaceType(i)=1;
                        end
                    end
                end
            end
        end
    end
end