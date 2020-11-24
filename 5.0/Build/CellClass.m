classdef CellClass
   %% cell clas
    properties
     X  
     Type                  %% 1 --> ag 2 ---> mono
     Ext                   %% (array-structure)IDs of extenal Cells  
     Int                   %% (array-structure)IDs of internal Cells 
     Tris                  %% (cell-structure) of matrices [ntris 3] with Triangles defining cell surfaces [surface-center v1 v1].
     cTet                  %% (cell-structure) of matrices [ntet 4]  with tetrahedrons connected to node i.         
     cTetID                %% (cell-structure) of arrays [ntet 1]    with IDs of the connected tetrahedrons. 
     cTri                  %% (cell-structure) of matrices [nTri 3]  with Tringels connected to node i.         
     cTriID                %% (cell-structure) of arrays [nTri 1]   with IDs of the connected tringles.
     cNodes                %% (cell-structure) of arrays [ntet 1]   with IDs of the connected nodes. 
     Cv                    %% (cell-structure) of matrices [nBars 2].
     CvID                  %% (cell-structure) of arrays [nBars 1] with IDs of cellular-bars.
     SurfsCenters          %% (matrix-structure) [nSurfaces 3] coordinate of surface centers.
%      CellSurfsCentersID  %% (cell-structure) of arrays  [nCellularSurfaces 1] IDs of surface-centers.
     nTotalTris            %% (scalar) total number of surfaces-triangles 
     n                     %% (scalar) number of cells 
     Vol                   %% (array-structure) [nCells 1] the volume of cells 
     Vol0                  %% (array-structure) [nCells 1] the volume of cells 
     SArea                 %% (array-structure) [nCells 1] the Surface Area of cells 
     SArea0                %% (array-structure) [nCells 1] the initial the Surface Area of cells 
     Surfaces              %% (cell-structure) {nCells 1} of cell-structure of arrays 
%                             .nSurfaces scaler  number of Surfaces; 
%                             .SurfaceCentersID  array [nCellularSurfaces 1] IDs of surface-centers. 
%                             .SurfaceVertices (cell-structure) of arrays  [nSurf-vertices 1] IDs of surface-vertices. 
     SAreaTri                 %% (array-structure) [nCells 1] 
     SAreaTri0                %% (array-structure) [nCells 1] 
     SAreaFace            %%  (array-structure) [nCells 1] 
     SAreaFace0           %%  (array-structure) [nCells 1] 
     AssembleAll            %%  logical 
     AssembleNodes          %%  (array-structure) list f node which correspond to the cell-center of the cells to be assmbled 
                            % it has no effect if (AssembleAll = true)
    Edges                        
    end
   methods
      function Cell = CellClass(X,nC,xInternal,xExternal)
         if nargin > 0 
            Cell.X = X;
            Cell.Ext=xExternal;
            Cell.Int=xInternal;
            Cell.Tris=cell(nC,1);
            Cell.Type=zeros(nC,1);
            Cell.cTet=cell(nC,1);
            Cell.cTetID=cell(nC,1);
            Cell.cTri=cell(nC,1);
            Cell.cTriID=cell(nC,1);
            Cell.cNodes=cell(nC,1);
            Cell.Cv=cell(nC,1);
            Cell.CvID=cell(nC,1);
            Cell.n=length(xInternal);
            Cell.Vol=zeros(Cell.n,1);
            Cell.Vol0=zeros(Cell.n,1);
            Cell.SArea=zeros(Cell.n,1);
            Cell.SArea0=zeros(Cell.n,1);
            Cell.SAreaTri=cell(Cell.n,1);
            Cell.SAreaTri0=cell(Cell.n,1);
            Cell.SAreaFace=cell(Cell.n,1);
            Cell.SAreaFace0=cell(Cell.n,1);
            Cell.Surfaces=cell(nC,1);
            Cell.AssembleAll=true;
            Cell.AssembleNodes=[];
            Cell.Edges=cell(nC,1);

         end
      end
      
      %% Ablate cells
      % We assume there will only be one ablation. Thus, we could remove
      % the IDs from the middle.
      % TODO: Check if Cell.Int need to be in order in consecutive numbers
      function cell = AblateCells(obj, cellsToRemove)
          obj.Int(cellsToRemove) = [];
          obj.Ext = [obj.Ext cellsToRemove];
          
          %Consequences:
          obj.Type(cellsToRemove) = [];
          %cell structures
          obj.Tris(cellsToRemove) = [];
          obj.cTet(cellsToRemove) = [];
          obj.cTri(cellsToRemove) = [];
          obj.cTetID(cellsToRemove) = [];
          obj.cTriID(cellsToRemove) = [];
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
          obj.SAreaTri0(cellsToRemove) = [];
          obj.SAreaFace(cellsToRemove) = [];
          obj.SAreaFace0(cellsToRemove) = [];
          obj.Surfaces(cellsToRemove) = [];
          obj.Edges(cellsToRemove) = [];
          
          obj.n = obj.n - length(cellsToRemove);
          obj.nTotalTris = sum(cellfun(@(X) size(X,1), obj.Tris));
          
          %Return
          cell = obj;
      end
      
      function cell = computeEdgeLengths(obj, Y)
          
          obj.EdgeLengths = [];
      end
   end
end