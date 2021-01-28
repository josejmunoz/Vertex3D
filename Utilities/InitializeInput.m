function [CellInput]=InitializeInput(Cell,Set)
% Initialize Lmabda-reduction-factor for each cell-inerface (check KgSurfaceCellBasedAdhesionParallel.m)
% Initialize propulsion force for each cell (check gPropulsion.m)

if Set.SurfaceType==4
    aux1=Set.LambdaS1CellFactor;
    aux2=Set.LambdaS2CellFactor;
    aux3=Set.LambdaS3CellFactor;
    aux4=Set.LambdaS4CellFactor; 
    CellInput.LambdaS1Factor=ones(Cell.n,1);
    CellInput.LambdaS2Factor=ones(Cell.n,1);
    CellInput.LambdaS3Factor=ones(Cell.n,1);
    CellInput.LambdaS4Factor=ones(Cell.n,1);
    if ~isempty(aux1)
        for i=1:size(aux1,1)
            CellInput.LambdaS1Factor(aux1(i,1))=aux1(i,2);
        end
    end
    
    if ~isempty(aux2)
        for i=1:size(aux2,1)
            CellInput.LambdaS2Factor(aux2(i,1))=aux2(i,2);
        end
    end
    
    if ~isempty(aux3)
        for i=1:size(aux3,1)
            CellInput.LambdaS3Factor(aux3(i,1))=aux3(i,2);
        end
    end
    
    if ~isempty(aux4)
        for i=1:size(aux4,1)
            CellInput.LambdaS4Factor(aux4(i,1))=aux4(i,2);
        end
    end   
else 
    CellInput.LambdaS1Factor=[];
    CellInput.LambdaS2Factor=[];
    CellInput.LambdaS3Factor=[];
    CellInput.LambdaS4Factor=[];
end 



if Set.Propulsion
    CellInput.Propulsion=zeros(Cell.n,3);
    CellInput.Propulsion(:,1:2)=(rand(Cell.n,2)-0.5)./50;
else 
    CellInput.Propulsion=[];
end 


end 
