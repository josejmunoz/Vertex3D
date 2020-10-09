function [g,Cell]=gVolume(Cell,Y,Set)
%% Set parameters
ncell=Cell.n;
%% Initialize
dimg=Set.NumTotalV*3;
g=zeros(dimg,1); % Local cell residual



%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);

%% Loop over Cells 
%     % Analytical residual g and Jacobian K
for i=1:ncell
    lambdaV=Set.lambdaV;
    fact=lambdaV*(Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i)^2;
    ge=zeros(dimg,1); % Local cell residual
    
%     YY=unique(Cell.Tris{i}(:,[1 2]));
%     YYY=unique(Cell.Tris{i}(Cell.Tris{i}(:,3)<0,3));
%     YYY=abs(YYY);
%     CC=unique(Cell.Tris{i}(Cell.Tris{i}(:,3)>0,3));
%     A=length(YY)+length(YYY)+length(CC);
%     YY=sum(Y.DataRow(YY,:),1);
%     YYY=sum(Y.DataRow(YYY,:),1);
%     CC=sum(Cell.SurfsCenters.DataRow(CC,:),1);
%     YY=(YY+YYY+CC)./A;
%     
    
    
    
    % Loop over Cell-face-triangles
    Tris=Cell.Tris{i};
    for t=1:size(Tris,1)
        nY=Tris(t,:);
        Y1=Y.DataRow(nY(1),:);
        Y2=Y.DataRow(nY(2),:);
        if nY(3)<0
            nY(3)=abs(nY(3));
            Y3=Y.DataRow(nY(3),:);
        else 
            Y3=Cell.SurfsCenters.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
        end 
        Y3=Y3;
        [gs]=gKDet(Y1,Y2,Y3);
        ge=AssemblegTriangleVol(ge,gs,nY);
    end 
    g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
end


end
%%
%%
function [gs]=gKDet(Y1,Y2,Y3)
% Returns residual and  Jacobian of det(Y)=y1'*cross(y2,y3)
% gs=[der_y1 det(Y) der_y2 det(Y) der_y3 det(Y)]
gs=[cross(Y2,Y3)'; % der_Y1 (det(Y1,Y2,Y3)) 
    cross(Y3,Y1)';
    cross(Y1,Y2)'];
end
%%
function   ge=AssemblegTriangleVol(ge,gt,nY)
% Assembles volume residual of a triangle of vertices (9 components)
dim=3;
for I=1:length(nY) % loop on 3 vertices of triangle
    if nY(I)>0
         idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
         idofl=(I-1)*dim+1:I*dim;
         ge(idofg)=ge(idofg)+gt(idofl);
    end
end
end
%%
%%



