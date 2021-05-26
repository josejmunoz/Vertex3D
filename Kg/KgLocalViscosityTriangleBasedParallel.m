function [g,K,Cell,EnergyS]=KgLocalViscosityTriangleBasedParallel(Cell,Y,Set)
% Local viscous effect based on the Area of Triangles
% Potential:  -> Set.LocalViscosityOption=1 -> W=(1/2) nu_Local/dt sum( ((At-Atn)/Atn)^2 )
%             -> Set.LocalViscosityOption=2 -> W=(1/2) nu_Local/dt sum( ((At-Atn))^2 )
%                     where At: is the Area of triangle at t_(n+1)
%                           Atn: is the Area of triangle t_(n)
%                           dt: time step


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

K=sparse(zeros(dimg)); % Also used in sparse

EnergyS=0;

CellAux.AssembleAll=Cell.AssembleAll;
CellAux.Int=Cell.Int;
CellAux.AssembleNodes=Cell.AssembleNodes;
CellAux.FaceCentres=Cell.FaceCentres;
CellAux.Tris=Cell.Tris;
CellAux.SAreaTri=Cell.SAreaTri;
CellAux.SAreaTrin=Cell.SAreaTrin;
CellAux.RemodelledVertices=Cell.RemodelledVertices;
CellAux.DebrisCells = Cell.DebrisCells;

%% Loop over Cells
%     % Analytical residual g and Jacobian K
parfor i=1:ncell
    if CellAux.DebrisCells(i)
        continue;
    end 
    if ~CellAux.AssembleAll
        if ~ismember(CellAux.Int(i),CellAux.AssembleNodes)
            continue
        end
    end
    lambdaS=Set.nu_Local_TriangleBased/Set.dt;
    % Loop over Cell-face-triangles
    Tris=CellAux.Tris{i};
    for t=1:size(Tris,1)
        ge=zeros(dimg,1); % Local cell residual
        if  Set.LocalViscosityOption==1
            fact=lambdaS *  (CellAux.SAreaTri{i}(t)-CellAux.SAreaTrin{i}(t)) / CellAux.SAreaTrin{i}(t)^2   ;
        elseif Set.LocalViscosityOption==2
            fact=lambdaS *  (CellAux.SAreaTri{i}(t)-CellAux.SAreaTrin{i}(t))    ;
        end
        nY=Tris(t,:);
        Y1=Y.DataRow(nY(1),:);
        Y2=Y.DataRow(nY(2),:);
        if nY(3)<0
            nY(3)=abs(nY(3));
            Y3=Y.DataRow(nY(3),:);
        else
            Y3=CellAux.FaceCentres.DataRow(nY(3),:);
            nY(3)=nY(3)+Set.NumMainV;
        end
        if ~CellAux.AssembleAll && ~any(ismember(nY,CellAux.RemodelledVertices))
            continue
        end
        [gs,Ks,Kss]=gKSArea(Y1,Y2,Y3);
        Ke=fact*(Ks+Kss);
        ge=Assembleg(ge,gs,nY);
        [si{i},sj{i},sv{i},sk{i}]= AssembleKSparse(Ke,nY,si{i},sj{i},sv{i},sk{i});
        g=g+ge*fact; 
        if  Set.LocalViscosityOption==1
            K=K+sparse((ge)*(ge')*lambdaS/(CellAux.SAreaTrin{i}(t)^2));
        elseif Set.LocalViscosityOption==2
            K=K+sparse((ge)*(ge')*lambdaS);
        end
        EnergyS=EnergyS+0;
    end
end

[si,sj,sv]=sReduction(si,sj,sv,sk,Cell);


K=sparse(si,sj,sv,dimg,dimg)+K;

end
%%
