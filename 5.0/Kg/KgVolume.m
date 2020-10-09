function [g,K,Cell,EnergyV]=KgVolume(Cell,Y,Set)
% K(i,j)= derivative of g(i) wrt to x(j)



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
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells 
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
else
    K=zeros(dimg); % Also used in sparse
end

EnergyV=0;
%% Compute Volume
[Cell]=ComputeCellVolume(Cell,Y);


%% Loop over Cells 
%     % Analytical residual g and Jacobian K
for i=1:ncell
    lambdaV=Set.lambdaV;
    fact=lambdaV*(Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i)^2;
    ge=zeros(dimg,1); % Local cell residual
    
    YY=unique(Cell.Tris{i}(:,[1 2]));
    YYY=unique(Cell.Tris{i}(Cell.Tris{i}(:,3)<0,3));
    YYY=abs(YYY);
    CC=unique(Cell.Tris{i}(Cell.Tris{i}(:,3)>0,3));
    A=length(YY)+length(YYY)+length(CC);
    YY=sum(Y.DataRow(YY,:),1);
    YYY=sum(Y.DataRow(YYY,:),1);
    CC=sum(Cell.SurfsCenters.DataRow(CC,:),1);
    YY=(YY+YYY+CC)./A;
    
    
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
        [gs,Ks]=gKDet(Y1,Y2,Y3);
        ge=AssemblegTriangleVol(ge,gs,nY);
        if Set.Sparse
            [si,sj,sv,sk]= AssembleKTriangleVolSparse(Ks*fact/6,nY,si,sj,sv,sk);
        else
            K= AssembleKTriangleVol(K,Ks*fact/6,nY);
        end
    end 
 
    g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
    if Set.Sparse
        K=K+lambdaV*sparse((ge)*(ge'))/6/6/Cell.Vol0(i)^2;
    else
          K=K+lambdaV*(ge)*(ge')/6/6/Cell.Vol0(i)^2;
    end
    EnergyV=EnergyV+ lambdaV/2 *((Cell.Vol(i)-Cell.Vol0(i))/Cell.Vol0(i))^2;

end

if Set.Sparse
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end
%%
function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];

end
%%
function [gs,Ks]=gKDet(Y1,Y2,Y3)
% Returns residual and  Jacobian of det(Y)=y1'*cross(y2,y3)
% gs=[der_y1 det(Y) der_y2 det(Y) der_y3 det(Y)]
% Ks=[der_y1y1 det(Y) der_y1y2 det(Y) der_y1y3 det(Y)
%     der_y2y1 det(Y) der_y2y2 det(Y) der_y2y3 det(Y)
%     der_y3y1 det(Y) der_y3y2 det(Y) der_y3y3 det(Y)]
dim=length(Y1);
gs=[cross(Y2,Y3)'; % der_Y1 (det(Y1,Y2,Y3)) 
    cross(Y3,Y1)';
    cross(Y1,Y2)'];
Ks=[ zeros(dim) -Cross(Y3)   Cross(Y2) % g associated to der wrt vertex 1
    Cross(Y3)   zeros(dim) -Cross(Y1)
    -Cross(Y2)   Cross(Y1)  zeros(dim)];
end
%%
function   ge=AssemblegTriangleVol(ge,gt,nY)
% Assembles volume residual of a triangle of vertices (9 components)
dim=3;
for I=1:length(nY) % loop on 3 vertices of triangle
    if nY(I)>0
         idofg=(nY(I)-1)*dim+1:nY(I)*dim;% global dof
         idofl=(I-1)*dim+1:I*dim;
         ge(idofg)=ge(idofg)+gt(idofl);
    end
end
end
%%
function Ke= AssembleKTriangleVol(Ke,Kt,nY)
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
function [si,sj,sv,sk]= AssembleKTriangleVolSparse(Kt,nY,si,sj,sv,sk)
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



