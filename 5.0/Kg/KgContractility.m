function [g,K,Cell,energy] = KgContractility(Cell,Y,Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

%% Set parameters
totalCells = Cell.n;
C = Set.cContractility;

%% Initialize
dimensionsG = Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells 
    sj=si;
    sv=si;
%     K=sparse(zeros(dimg)); % Also used in sparse
else
    K=zeros(dimg); % Also used in sparse
end

energy = 0;

%% Calculate basic information
[Cell] = Cell.computeEdgeLengths();
Cell.EdgeLengths = Cell.EdgeLengths;
edgeLengths0 = Cell.EdgeLengths0 * C;

%% Calculate residual g
g = computeGContractility(Cell, C);

%% Calculate Jacobian
K = computeKContractility(Cell, C);

%% Calculate energy
energy = energy + computeEnergyContractility(edgeLength, edgeLength0, C);

end

function [kContractility] = computeKContractility(Cell, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

dim = 3;

kContractility = (1 - (edgeLength0)/edgeLength) * eye(dim) + (edgeLength0/edgeLength^2) * (1/edgeLength) * cross(y_1 - y_2, y_1 - y_2);

end

function [gContractility] = computeGContractility(Cell, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

gContractility = (edgeLength - edgeLength0) * (y_1 - y_2) / edgeLength;

end

function [energyConctratility] = computeEnergyContractility(edgeLength, edgeLength0, C)
    energyConctratility = 0.5 * (edgeLength - (edgeLength0 * C))^2;
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
function Ke= AssembleKTriangleSArea(Ke,Kt,nY)
% Assembles volume Jacobian of a triangle of vertices (9x9 components)
dim=3;


for I=1:length(nY) % loop on 3 vertices of triangle
    for J=1:length(nY)
        if nY(I)>0
            idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
            idofl=(I-1)*dim+1:I*dim;
            jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
            jdofl=(J-1)*dim+1:J*dim;
            Ke(idofg,jdofg)=Ke(idofg,jdofg)+Kt(idofl,jdofl);
        end
    end 
end
end


%%
function [si,sj,sv,sk]= AssembleKTriangleSAreaSparse(Kt,nY,si,sj,sv,sk)
% Assembles volume Jacobian of a triangle of vertices (9x9 components)
dim=3;

for I=1:length(nY) % loop on 3 vertices of triangle
    for J=1:length(nY)
        if nY(I)>0
            idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
            idofl=(I-1)*dim+1:I*dim;
            jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
            jdofl=(J-1)*dim+1:J*dim;
            for d=1:dim
                si(sk+1:sk+dim)=idofg;
                sj(sk+1:sk+dim)=jdofg(d);
                sv(sk+1:sk+dim)=Kt(idofl,jdofl(d));
                sk=sk+dim;
            end
        end
    end 
end

end
%%