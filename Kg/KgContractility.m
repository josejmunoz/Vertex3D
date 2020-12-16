function [g,K,Cell,energy] = KgContractility(Cell,Y,Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

%% Set parameters
C = Set.cContractility;
Set.Sparse = true;

%% Initialize
dimensionsG = Set.NumTotalV*3;
g=zeros(dimensionsG,1); % Local cell residual

if Set.Sparse
    sk=0;
    si=zeros((dimensionsG*3)^2,1); % Each vertex is shared by at least 3 cells 
    sj=si;
    sv=si;
    K=sparse(zeros(dimensionsG)); % Also used in sparse
else
    K=zeros(dimensionsG); % Also used in sparse
end

energy = 0;

%% Calculate basic information
[Cell] = Cell.computeEdgeLengths(Y);

for numCell = 1:Cell.n
    if ~Cell.GhostCells(numCell)
        continue;
    end
    
    edgeLengths = Cell.EdgeLengths{numCell};
    edgeLengths0 = Cell.EdgeLengths0{numCell};
    edgeLengths0(:, 1) = edgeLengths0(:, 1) * C;
    
    for numEdge = 1:size(edgeLengths, 3)
        
        y_1 = Y.DataRow(edgeLengths(numEdge, 2), :);
        y_2 = Y.DataRow(edgeLengths(numEdge, 3), :);
        
        l_i = edgeLengths(numEdge, 1);
        l_i0 = edgeLengths0(numEdge, 1);
        
        %% Calculate residual g
        g_current = computeGContractility(l_i0, l_i, y_1, y_2);
        g = assembleG(g, g_current, edgeLengths(numEdge, 2:3));
        %% Calculate Jacobian
        K_current = computeKContractility(l_i0, l_i, y_1, y_2);
        
        if Set.Sparse
            [si,sj,sv,sk] = assembleKSparse(K_current, edgeLengths(numEdge, 2:3), si, sj, sv, sk);
        else
            K = assembleK(K, K_current, edgeLengths(numEdge, 2:3));
        end
        
        %% Calculate energy
        energy = energy + computeEnergyContractility(l_i, l_i0);
    end
    
    if Set.Sparse
        K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
    end
end
end

function [kContractility] = computeKContractility(l_i0, l_i, y_1, y_2)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

dim = 3;

kContractility(1:3, 1:3) = (1 - (l_i0/l_i)) * eye(dim) + (l_i0/l_i^2) * (1/l_i) * (y_1 - y_2).*(y_1 - y_2)';
kContractility(1:3, 4:6) = -kContractility(1:3, 1:3);
kContractility(4:6, 1:3) = -kContractility(1:3, 1:3);
kContractility(4:6, 4:6) = kContractility(1:3, 1:3);

end

function [gContractility] = computeGContractility(l_i0, l_i, y_1, y_2)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

gContractility(1:3, 1) = (l_i - l_i0) * (y_1 - y_2) / l_i;
gContractility(4:6, 1) = -gContractility(1:3);

end

function [energyConctratility] = computeEnergyContractility(l_i, l_i0)
    energyConctratility = 0.5 * (l_i - (l_i0))^2;
end

%%
function   gAccumulated=assembleG(gAccumulated,gCurrent,nY)
% Assembles residual of a pair of vertices (6 components)
dim=3;
for I=1:length(nY) % loop on 2 vertices of an edge
    if nY(I)>0
         idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
         idofl=(I-1)*dim+1:I*dim;
         gAccumulated(idofg)=gAccumulated(idofg)+gCurrent(idofl);
    end
end
end
%%
function Kaccumulated= assembleK(Kaccumulated,Kcurrent,nY)
% Assembles volume Jacobian of a triangle of vertices (12x12 components)
dim=3;


for I=1:length(nY)
    for J=1:length(nY)
        if nY(I)>0
            idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
            idofl=(I-1)*dim+1:I*dim;
            jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
            jdofl=(J-1)*dim+1:J*dim;
            Kaccumulated(idofg,jdofg)=Kaccumulated(idofg,jdofg)+Kcurrent(idofl,jdofl);
        end
    end 
end
end


%%
function [si,sj,sv,sk]= assembleKSparse(Kcurrent,nY,si,sj,sv,sk)
% Assembles Jacobian of a pair of vertices (6x6 components)
dim=3;

for I=1:length(nY)
    for J=1:length(nY)
        if nY(I)>0
            idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
            idofl=(I-1)*dim+1:I*dim;
            jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
            jdofl=(J-1)*dim+1:J*dim;
            for d=1:dim
                si(sk+1:sk+dim)=idofg;
                sj(sk+1:sk+dim)=jdofg(d);
                sv(sk+1:sk+dim)=Kcurrent(idofl,jdofl(d));
                sk=sk+dim;
            end
        end
    end 
end

end
%%