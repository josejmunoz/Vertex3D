function [g,K,Cell,EnergyS]=KgSurfaceCellBasedAdhesion(Cell,Y,Faces,Set)
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
if Set.Sparse
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells 
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
else
    K=zeros(dimg); % Also used in sparse
end


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
    fact0=0;
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

        fact0=fact0+Lambda*Cell.SAreaFace{i}(f);
    end 
    fact=fact0/Cell.SArea0(i)^2;

    
    
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
            [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
            Ks=fact*Lambda*(Ks+Kss);
            gs=Lambda*gs;
            ge=AssemblegTriangleSArea(ge,gs,nY);
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKTriangleSAreaSparse(Ks,nY,si,sj,sv,sk);
            else
                K= AssembleKTriangleSArea(K,Ks,nY);
            end
        end

    end 
    
    g=g+ge*fact; % Volume contribution of each triangle is (y1-y2)'*J*(y2-y3)/2
    if Set.Sparse
        K=K+sparse((ge)*(ge')/(Cell.SArea0(i)^2));
    else
        K=K+(ge)*(ge')/(Cell.SArea0(i)^2);  %-(gee)*(gee')*(fact);
    end
    
    
    % Cell.SAreaFace0 need to be updated  after Remodling 
     EnergyS=EnergyS+ (1/2)*fact0*fact;

end

if Set.Sparse
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

