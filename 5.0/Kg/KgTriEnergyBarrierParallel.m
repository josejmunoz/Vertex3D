function [g,K,Cell,EnergyB]=KgTriEnergyBarrierParallel(Cell,Y,Set)
% K(i,j)= derivative of g(i) wrt to x(j)
% energy based on the local cell area  (triangles) W_s= sum_triangles ((At-At0)/At0)^2

%% Input 
Set.Sparse=true;
% Set.Sparse=false;

% Set.yRelaxtion Main
% Set.nodes
% Set.nv
% Set.lambdaV

%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse
    si=cell(ncell,1); %zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells 
    sj=si;
    sv=si;
    sk=si;
    K=sparse(zeros(dimg)); % Also used in sparse
else
    K=zeros(dimg); % Also used in sparse
end


for i=1:ncell
    si{i}=zeros(size(Cell.Edges{i}*3*3,1),1);
    sj{i}=si{i};
    sv{i}=si{i};
    sk{i}=0;
end



EnergyB=0;
%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);
Beta=Set.Beta;
lambdaB=Set.lambdaB;
NumMainV=Set.NumMainV;
CellAssembleAll=Cell.AssembleAll;
CellInt=Cell.Int;
CellAssembleNodes=Cell.AssembleNodes;
CellTris=Cell.Tris;
CellSAreaTri=Cell.SAreaTri;
CellSAreaTri0=Cell.SAreaTri0;
CellSurfsCenters=Cell.SurfsCenters;

%% Loop over Cells 
%     % Analytical residual g and Jacobian K
parfor i=1:ncell
    if ~CellAssembleAll
        if ~ismember(CellInt(i),CellAssembleNodes) 
           continue
        end 
    end 
    % Loop over Cell-face-triangles
    Tris=CellTris{i};
    for t=1:size(Tris,1)
%         Cell.SAreaTri0{i}(t)=1e-4;
       fact=-((lambdaB*Beta)/CellSAreaTri0{i}(t)) * exp(lambdaB*(1-Beta*CellSAreaTri{i}(t)/CellSAreaTri0{i}(t)));
       fact2=fact*-((lambdaB*Beta)/CellSAreaTri0{i}(t));
        nY=Tris(t,:);
        Y1=Y.DataRow(nY(1),:); %#ok<PFBNS>
        Y2=Y.DataRow(nY(2),:);
        Y3=CellSurfsCenters.DataRow(nY(3),:); %#ok<PFBNS>
        nY(3)=nY(3)+NumMainV;
        [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
        Ks=(gs)*(gs')*fact2+Ks*fact+Kss*fact; 
%         if Set.Sparse
            [si{i},sj{i},sv{i},sk{i}]= AssembleKTriangleSAreaSparse(Ks,nY,si{i},sj{i},sv{i},sk{i});
%         else
%             K= AssembleKTriangleSArea(K,Ks,nY);
%         end
        gaux=zeros(dimg,1);
        gs=fact.*gs;
        dim=3;
        for I=1:length(nY) % loop on 3 vertices of triangle
             idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
             idofl=(I-1)*dim+1:I*dim;
             gaux(idofg)=gaux(idofg)+gs(idofl);
        end
        g=g+gaux
        EnergyB=EnergyB+ exp(lambdaB*(1-Beta*CellSAreaTri{i}(t)/CellSAreaTri0{i}(t)));
    end 
end



[si,sj,sv]=sReduction(si,sj,sv,sk,Cell);

% if Set.Sparse
    K=sparse(si,sj,sv,dimg,dimg)+K;
% end

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

function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];
end


function [Si,Sj,Sv]=sReduction(si,sj,sv,sk,Cell)
Sk=0;
for i=1:Cell.n
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end
    Sk=Sk+sk{i};
end 
Si=zeros(Sk,1); % Each vertex affecting 6 nodes
Sj=Si;
Sv=Si;
Sk=0;
for i=1:length(si)
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end
    Si(Sk+1:Sk+sk{i})=si{i};
    Sj(Sk+1:Sk+sk{i})=sj{i};
    Sv(Sk+1:Sk+sk{i})=sv{i};
    Sk=Sk+sk{i};
end 
end 
