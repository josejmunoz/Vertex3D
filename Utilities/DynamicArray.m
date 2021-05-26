classdef DynamicArray 
   %% cell clas
    properties
     DataRow  % [n 2] the Ids of the cells sharing this (i) face 
     NotEmpty  % [n 2] the Ids of the nodes sharing this (i) face 
     EmptyList     % {n 1} the tets corresponding to this face 
     nE        % counters on ListEmpty list
     n        % number of faces
    end
   methods
      function Array= DynamicArray(S1,S2)
          Array.n=0;
          Array.DataRow=zeros(S1,S2);
          Array.NotEmpty=false(S1,1);
          Array.EmptyList=zeros(S1,1);
          Array.nE=0;
      end 
      %% Size
      function s=Size(obj)  
               s=obj.n-obj.nE;
      end 
      %% Stored data
      function SD=Data(obj)
          SD=obj.DataRow(obj.NotEmpty,:);
      end 
      
      %% total Data
      function TD=DataOrdered(obj)
          TD=obj.DataRow(1:obj.n,:);
      end 
      %% Assign Data
      function [obj]=Modify(obj,D)
         obj.DataRow(1:size(D,1),:)=D;
      end 
      %% Add Data
      function [obj,list] = Add(obj,D)
          if size(D,2)~=size(obj.DataRow,2)
              error('error using DynamicArray.Add')
          end 
          list=zeros(size(D,1),1);
         for i=1:size(D,1)
             if obj.nE>0
                 obj.DataRow(obj.EmptyList(obj.nE),:)=D(i,:);
                 list(i)=obj.EmptyList(obj.nE);
                 obj.NotEmpty(obj.EmptyList(obj.nE))=true;
                 obj.EmptyList(obj.nE)=0;
                 obj.nE=obj.nE-1;
             else
                 obj.DataRow( obj.n+1 : obj.n+size(D(i:end,:),1) ,:)=D(i:end,:);
                 list(i:end)=obj.n+1 : obj.n+size(D(i:end,:),1);
                 obj.NotEmpty( obj.n+1 : obj.n+size(D(i:end,:),1) )=true;
                 obj.EmptyList( obj.n+1 : obj.n+size(D(i:end,:),1) ) = 0;
                 obj.n=obj.n+size(D(i:end,:),1);
                 break
             end 
         end
      end 
      %% Remove Data 
      function [obj] = Remove(obj,V)
               obj.DataRow(V,:)=0;
               obj.NotEmpty(V)=false;
               obj.EmptyList(obj.nE+1:obj.nE+length(V))=V;
               obj.nE=obj.nE+length(V);
      end
      
      
   end
end