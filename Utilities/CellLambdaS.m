function [Set]=CellLambdaS(Celln,Set)
if Set.SurfaceType~=4
    return
end 
    
aux1=Set.LambdaS1CellFactor;
aux2=Set.LambdaS2CellFactor;
aux3=Set.LambdaS3CellFactor; 

Set.LambdaS1CellFactor=ones(Celln,1);
Set.LambdaS2CellFactor=ones(Celln,1);
Set.LambdaS3CellFactor=ones(Celln,1);

if ~isempty(aux1)
    for i=1:size(aux1,1)
        Set.LambdaS1CellFactor(aux1(i,1))=aux1(i,2);
    end 
end

if ~isempty(aux2)
    for i=1:size(aux2,1)
        Set.LambdaS2CellFactor(aux2(i,1))=aux2(i,2);
    end 
end

if ~isempty(aux3)
    for i=1:size(aux3,1)
        Set.LambdaS3CellFactor(aux3(i,1))=aux3(i,2);
    end 
end 

end 