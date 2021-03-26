function [g,K,Cell,EnergyB]=KgTriEnergyBarrierParallel(Cell,Y,Set)
% The residual g and Jacobian K of  Energy Barrier
% Energy  WBexp = exp( Set.lambdaB*  ( 1 - Set.Beta*At/Set.BarrierTri0 )  );


%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
si=cell(ncell,1); %zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
sj=si;
sv=si;
sk=si;
for i=1:ncell
    si{i}=zeros(size(Cell.Edges{i},1)*3*3*10,1);
    sj{i}=si{i};
    sv{i}=si{i};
    sk{i}=0;
end

K=sparse(zeros(dimg)); 
EnergyB=0;
Beta=Set.Beta;
lambdaB=Set.lambdaB;
CellAssembleAll=Cell.AssembleAll;
CellInt=Cell.Int;
CellAssembleNodes=Cell.AssembleNodes;
CellTris=Cell.Tris;
CellSAreaTri=Cell.SAreaTri;
CellFaceCentres=Cell.FaceCentres;
CellRemodelledVertices=Cell.RemodelledVertices;
CellDebrisCells = Cell.DebrisCells;

%% Loop over Cells
% Analytical residual g and Jacobian K
parfor i=1:ncell
    if CellDebrisCells(i)
        continue;
    end 
    if ~CellAssembleAll
        if ~ismember(CellInt(i),CellAssembleNodes)
            continue
        end
    end
    % Loop over Cell-face-triangles
    Tris=CellTris{i};
    for t=1:size(Tris,1)
        fact=-((lambdaB*Beta)/Set.BarrierTri0) * exp(lambdaB*(1-Beta*CellSAreaTri{i}(t)/Set.BarrierTri0));
        fact2=fact*-((lambdaB*Beta)/Set.BarrierTri0);
        nY=Tris(t,:);
        Y1=Y.DataRow(nY(1),:); %#ok<PFBNS>
        Y2=Y.DataRow(nY(2),:);
        if nY(3)<0
            nY(3)=abs(nY(3));
            Y3=Y.DataRow(nY(3),:);
        else
            Y3=CellFaceCentres.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
        end
        if ~CellAssembleAll && ~any(ismember(nY,CellRemodelledVertices))
            continue
        end
        [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
        Ks=(gs)*(gs')*fact2+Ks*fact+Kss*fact;
        [si{i},sj{i},sv{i},sk{i}]= AssembleKSparse(Ks,nY,si{i},sj{i},sv{i},sk{i});
        gaux=zeros(dimg,1);
        gs=fact.*gs;
        dim=3;
        for I=1:length(nY) % loop on 3 vertices of triangle
            idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
            idofl=(I-1)*dim+1:I*dim;
            gaux(idofg)=gaux(idofg)+gs(idofl);
        end
        g=g+gaux
        EnergyB=EnergyB+ exp(lambdaB*(1-Beta*CellSAreaTri{i}(t)/Set.BarrierTri0));
    end
end

[si,sj,sv]=sReduction(si,sj,sv,sk,Cell);

K=sparse(si,sj,sv,dimg,dimg)+K;

end
%%

%%
function [gs,Ks,Kss]=gKSArea(y1,y2,y3)
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


Kss=-(2/norm(q)).*(gs)*(gs');
Ks=fact.*[Q1'*Q1               KK(1,2,3,y1,y2,y3) KK(1,3,2,y1,y2,y3);
    KK(2,1,3,y1,y2,y3)  Q2'*Q2              KK(2,3,1,y1,y2,y3);
    KK(3,1,2,y1,y2,y3)  KK(3,2,1,y1,y2,y3) Q3'*Q3            ];
end
function [KIJ]=KK(i,j,k,y1,y2,y3)
Y=[y1;y2;y3];
%KIJ= (Yk-Yj)*(Yi-Yk)+Cross()
KIJ= (Cross(Y(j,:))-Cross(Y(k,:)))*(Cross(Y(i,:))-Cross(Y(k,:)))+...
    Cross(Cross(Y(j,:))*Y(i,:)')-Cross(Cross(Y(j,:))*Y(k,:)')-Cross(Cross(Y(k,:))*Y(i,:)');

end

