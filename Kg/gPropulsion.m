function [gp]=gPropulsion(Cell,Y,Set,CellInput,XgSub)




dimg=Set.NumTotalV*3;

gp=zeros(dimg,1); % Local cell residual




%% Loop over Cells 
%     % Analytical residual g and Jacobian K
for i=1:Cell.n
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end 
    for f=1:Cell.Surfaces{i}.nSurfaces
        if ismember(Cell.cNodes{i}(f),XgSub)
%             VSC=[Cell.Surfaces{i}.SurfaceVertices{f};...
%                  Cell.Surfaces{1}.SurfaceCentersID+Y.n]; 
             VSC=[Cell.Surfaces{i}.SurfaceVertices{f}];
            Dofs=3.*(kron(VSC,[1; 1; 1])-1)+kron(ones(length(VSC),1),[1; 2; 3]);
            PForcs=kron(ones(length(VSC),1),CellInput.Propulsion(i,:)');
            gp(Dofs)=PForcs;            
        end 
    end 
    
    
    
end 




end 