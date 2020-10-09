function [g,EnergyS]=gSurfaceCellBasedAdhesion(Cell,Y,Faces,Set)
% K(i,j)= derivative of g(i) wrt to x(j)
% energy based on the total cell area W_s= sum_cell ((As-As0)/As0)^2

%% Input 
Set.Sparse=true;
%  Set.Sparse=false;
% 
% Set.yRelaxtion Main
% Set.nodes
% Set.nv
% Set.lambdaV

%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual



EnergyS=0;

%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);

%% Loop over Cells 
%     % Analytical residual g and Jacobian K
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end 
    ge=zeros(dimg,1); % Local cell residual

    % first loop to commpute fact 
     fact=0;

    for f=1:Cell.Surfaces{i}.nSurfaces        
        if Faces.InterfaceType(Cell.Surfaces{i}.SurfaceCentersID(f))==0
            % External
            Lambda=Set.lambdaS1*Set.LambdaS1CellFactor(i);
        elseif  Faces.InterfaceType(Cell.Surfaces{i}.SurfaceCentersID(f))==1
            % Cell-Cell
            Lambda=Set.lambdaS2*Set.LambdaS2CellFactor(i);
        elseif Faces.InterfaceType(Cell.Surfaces{i}.SurfaceCentersID(f))==2
            % Cell-substrate
            Lambda=Set.lambdaS3*Set.LambdaS3CellFactor(i); 
        end
        % Loop over Cell-face-triangles
        Tris=Cell.Surfaces{i}.Tris{f};
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
            [gs]=gKSArea(Y1,Y2,Y3);
            gs=Lambda*gs;
            ge=AssemblegTriangleSArea(ge,gs,nY);
        end
        fact=fact+Lambda*Cell.SAreaFace{i}(f);
    end 
    fact=fact/Cell.SArea0(i)^2;

    g=g+ge*fact; % Volume contribution of each triangle is (y1-y2)'*J*(y2-y3)/2

%     EnergyS=EnergyS+ 0/2 *((Cell.SAreaFace{i}(f)-Cell.SAreaFace0{i}(f)) / Cell.SAreaFace0{i}(f))^2;

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

function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];

end

