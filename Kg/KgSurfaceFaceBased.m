function [g,K,Cell,EnergyS]=KgSurfaceFaceBased(Cell,Y,Set)
% The residual g and Jacobian K of Surface Energy
% Energy based on the area of faces W_s= lambdaS* sum_cell ((Af-Af0)/A0)^2

%% Input
Set.Sparse=true;


%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse && nargout>1
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
elseif nargout>1
    K=zeros(dimg); % Also used in sparse
end


EnergyS=0;

%% Loop over Cells
%     % Analytical residual g and Jacobian K
for i=1:ncell
    if Cell.DebrisCells(i)
        continue;
    end    
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
            Y1=Y.DataRow(nY(1),:);
            Y2=Y.DataRow(nY(2),:);
            Y3=Cell.FaceCentres.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
            [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
            ge=Assembleg(ge,gs,nY);
            if nargout>1
                Ke=Ks*fact+Kss*fact;
                if Set.Sparse
                    [si,sj,sv,sk]= AssembleKSparse(Ke,nY,si,sj,sv,sk);
                else
                    K= AssembleK(K,Ke,nY);
                end
            end
        end
        g=g+ge*fact; 
        if nargout>1
            if Set.Sparse
                K=K+sparse((ge)*(ge')*lambdaS/(Cell.SAreaFace0{i}(f)^2));
            else
                K=K+(ge)*(ge')*lambdaS/(Cell.SAreaFace0{i}(f)^2); 
            end
            EnergyS=EnergyS+ lambdaS/2 *((Cell.SAreaFace{i}(f)-Cell.SAreaFace0{i}(f)) / Cell.SAreaFace0{i}(f))^2;
        end
    end
end

if Set.Sparse &&  nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
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

