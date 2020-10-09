classdef FacesClass
   %% cell clas
    properties
     Energy
     EnergyTri
     Vertices
     Area
     Area0
     AreaTri
     AreaTri0
     Nodes
     InterfaceType
     NotEmpty  % [n 2] the Ids of the nodes sharing this (i) face 
     EmptyList     % {n 1} the tets corresponding to this face 
     V3
     V4
     nE        % counters on ListEmpty list
     n        % number of faces
    end
   methods
      function Array= FacesClass(S1)
          Array.Vertices=cell(S1,1);
          Array.AreaTri=cell(S1,1);
          Array.AreaTri0=cell(S1,1);
          Array.Area=zeros(S1,1);
          Array.Area0=zeros(S1,1);         
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
      %% Size
      function s=Size(obj)  
               s=obj.n-obj.nE;
      end 
      %% Actual data
      %% Add Data
      function [obj,list] = Add(obj,Nodes,Vertices,Y,SurfsCenters)
             if obj.nE>0
                 obj.Vertices{obj.EmptyList(obj.nE)}=Vertices; 
                 obj.Nodes(obj.EmptyList(obj.nE),:)=Nodes;
                 obj.Energy(obj.EmptyList(obj.nE),:)=0;
                 if length(Vertices)==3
                     obj.V3(obj.EmptyList(obj.nE))=true;
                     obj.V4(obj.EmptyList(obj.nE))=false;
                     % Compute Area
                     YTri=Y(Vertices,:);
                     obj.AreaTri{obj.EmptyList(obj.nE)}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                     obj.AreaTri0{obj.EmptyList(obj.nE)}=1e-3;
                      obj.EnergyTri{obj.EmptyList(obj.nE)}=0;
                     obj.Area(obj.EmptyList(obj.nE))=sum(obj.AreaTri{obj.EmptyList(obj.nE)});
                     obj.Area0(obj.EmptyList(obj.nE))=obj.Area(obj.EmptyList(obj.nE));
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
                     obj.AreaTri0{obj.EmptyList(obj.nE)}=ones(length(Vertices),1)*1e-3;
                     obj.EnergyTri{obj.EmptyList(obj.nE)}=zeros(length(Vertices),1);
                     for j=1:length(Vertices)-1
                            YTri=[Y(Vertices([j j+1]),:); SurfsCenters(obj.EmptyList(obj.nE),:)];
                            obj.AreaTri{obj.EmptyList(obj.nE)}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
    %                         obj.AreaTri0{obj.EmptyList(obj.nE)}(j)=obj.AreaTri{obj.EmptyList(obj.nE)}(j);
                     end
                     YTri=[Y(Vertices([j+1 1]),:); SurfsCenters(obj.EmptyList(obj.nE),:)];
                     obj.AreaTri{obj.EmptyList(obj.nE)}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));

                     obj.Area(obj.EmptyList(obj.nE))=sum(obj.AreaTri{obj.EmptyList(obj.nE)});
                     obj.Area0(obj.EmptyList(obj.nE))=obj.Area(obj.EmptyList(obj.nE));

                     % obj.AreaTri0{obj.EmptyList(obj.nE)}(j+1)=obj.AreaTri{obj.EmptyList(obj.nE)}(j+1);

                     list=obj.EmptyList(obj.nE);
                     obj.NotEmpty(obj.EmptyList(obj.nE))=true;
                     obj.EmptyList(obj.nE)=nan;
                     obj.nE=obj.nE-1;
                 end                  

             else
                 obj.Vertices{obj.n+1}=Vertices;
                 obj.Nodes(obj.n+1,:)=Nodes;
                 obj.Energy(obj.n+1)=0;
                  if length(Vertices)==3
                     obj.V3(obj.n+1)=true;
                     obj.V4(obj.n+1)=false;
                     % Compute Area
                     YTri=Y(Vertices,:);
                     obj.AreaTri{obj.n+1}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                     obj.AreaTri0{obj.n+1}=1e-3;
                     obj.EnergyTri{obj.n+1}=0;
                     obj.Area(obj.n+1)=sum(obj.AreaTri{obj.n+1});
                     obj.Area0(obj.n+1)=obj.Area(obj.n+1);
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
                     obj.AreaTri0{obj.n+1}=ones(length(Vertices),1)*1e-3;
                     obj.EnergyTri{obj.n+1}=zeros(length(Vertices),1);
                     for j=1:length(Vertices)-1
                            YTri=[Y(Vertices([j j+1]),:); SurfsCenters(obj.n+1,:)];
                            obj.AreaTri{obj.n+1}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
    %                            obj.AreaTri0{obj.n+1}(j)=obj.AreaTri{obj.n+1}(j);
                     end
                     YTri=[Y(Vertices([j+1 1]),:); SurfsCenters(obj.n+1,:)];
                     obj.AreaTri{obj.n+1}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));

                     obj.Area(obj.n+1)=sum(obj.AreaTri{obj.n+1});
                     obj.Area0(obj.n+1)=obj.Area(obj.n+1);

    %                  obj.AreaTri0{obj.n+1}(j+1)=obj.AreaTri{obj.n+1}(j+1);
                     list=obj.n+1;
                     obj.NotEmpty(obj.n+1)=true;
                     obj.n=obj.n+1;
                  end

             end 
      end 
      %% Remove Data 
      function [obj] = Remove(obj,V)
          for i=1:length(V)
              obj.Vertices{V(i)}=[];
          end
          obj.Nodes(V,:)=nan;
          obj.Energy(V)=nan;
          obj.Area(V)=nan;
          obj.Area0(V)=nan;
          obj.InterfaceType(V)=0;
          obj.AreaTri{V}=[];
          obj.AreaTri0{V}=[];
          obj.NotEmpty(V)=false;
          obj.V3(V)=false;
          obj.V4(V)=false;
          obj.EmptyList(obj.nE+1:obj.nE+length(V))=V;
          obj.nE=obj.nE+length(V);
      end
      
      %% Compute Face Area
        function [obj]=ComputeAreaTri(obj,Y,SurfsCenters)
           for i=1:obj.n
               if obj.NotEmpty(i)
                   if length(obj.Vertices{i})==3
                     YTri=Y(obj.Vertices{i},:);
                     obj.AreaTri{i}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                     obj.AreaTri0{i}=1e-3;
                     obj.Area(i)=sum(obj.AreaTri{i});
                   else 
                       obj.AreaTri{i}=zeros(length(obj.Vertices{i}),1);
                       obj.AreaTri0{i}=ones(length(obj.Vertices{i}),1)*1e-3;
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
                   obj.AreaTri0{i}=[];
               end
          end
      end 

      
      %% Compute Face Enrgy
      function [obj]=ComputeEnergy(obj,Set)
          for i=1:obj.n
            if obj.NotEmpty(i)
               if length(obj.Vertices{i})==3
                    obj.EnergyTri{i}=exp(Set.lambdaB*(1-Set.Beta*obj.AreaTri{i}/obj.AreaTri0{i}));
                    obj.Energy(i)=sum(obj.EnergyTri{i});
               else 
                    obj.EnergyTri{i}=zeros(size(obj.Vertices{i}));
                    for j=1:length(obj.Vertices{i})
                        obj.EnergyTri{i}(j)= exp(Set.lambdaB*(1-Set.Beta*obj.AreaTri{i}(j)/obj.AreaTri0{i}(j)));
                    end
                    obj.Energy(i)=sum(obj.EnergyTri{i});
               end 
            else 
                obj.EnergyTri{i}=[];
                obj.Energy(i)=0;
            end 
          end
      end 
      
      %% Check Interior faces  
      
      function [obj]=CheckInteriorFaces(obj,XgID,XgSub)
          if nargin==2
              for i=1:obj.n
                  if obj.NotEmpty(i)
                      if any(ismember(obj.Nodes(i,:),XgID))
                          % cell-cell face
                          obj.InterfaceType(i)=0;
                      else
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
                          % cell-subdtrate face
                          obj.InterfaceType(i)=2;
                      else 
                         % external face
                          obj.InterfaceType(i)=1;
                      end
                  end
              end
          end
      end

      
   end
end