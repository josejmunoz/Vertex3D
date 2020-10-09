function [g,Cell]=gSurfaceFaceBased(Cell,Y,Set)
% K(i,j)= derivative of g(i) wrt to x(j)
% energy based on the local cell area  (faces) W_s= sum_triangles ((Af-Af0)/Af0)^2

%% Input 
Set.Sparse=true;
Set.Sparse=false;

% Set.yRelaxtion Main
% Set.nodes
% Set.nv
% Set.lambdaV

%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual




%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);

%% Loop over Cells 
%     % Analytical residual g and Jacobian K
EnergyS=0;
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end 
    lambdaS=Set.lambdaS;
    
    for f=1:Cell.Surfaces{i}.nSurfaces
        ge=zeros(dimg,1); % Local cell residual
        fact=lambdaS *  (Cell.SAreaFace{i}(f)-Cell.SAreaFace0{i}(f)) / Cell.SAreaFace0{i}(f)^2   ;
        
        % Loop over Cell-face-triangles
        Tris=Cell.Surfaces{i}.Tris{f};
        for t=1:size(Tris,1)
            nY=Tris(t,:);
            Y1=Y(nY(1),:);
            Y2=Y(nY(2),:);
            Y3=Cell.SurfsCenters(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
            [gs]=gKSArea(Y1,Y2,Y3);
            ge=AssemblegTriangleSArea(ge,gs,nY);
        end
        g=g+ge*fact; % 
    end 
 end


end
%%

%%
function [gs]=gKSArea(y1,y2,y3)
% Returns residual and  Jacobian of A=||(y1-y2)x(y2-y3)||
% notation
% YI=Cross(yI)
% q =(y1-y3)x(y2-y3) = Y1y2 - Y1y3 - Y3y2
q= Cross(y2)*y1' - Cross(y2)*y3' + Cross(y1)*y3';
% Q_I = der_yI q =  Cross(yK) -Cross(yJ)
Q1=Cross(y2)-Cross(y3);
Q2=Cross(y3)-Cross(y1);
Q3=Cross(y1)-Cross(y2);
% KK_IJ = der_yJ QI
fact=1/(2*norm(q));
gs=fact.*[Q1'*q; % der_Y1 (det(Y1,Y2,Y3)) 
          Q2'*q;
          Q3'*q];
               
end



%%
function   ge=AssemblegTriangleSArea(ge,gt,nY)
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
%%

function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];

end

