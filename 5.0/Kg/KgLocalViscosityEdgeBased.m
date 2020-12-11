function [g,K,Cell,EnergyBar]=KgLocalViscosityEdgeBased(Cell,Y,Set)
% Local viscous effect based on the length of the edges between vertices
% Potential: -> Set.BendingAreaDependent=1 : W=(1/2) nu_Local/dt sum( ((L-Ln)/Ln)^2 )
%            -> Set.LocalViscosityOption=2 : W=(1/2) nu_Local/dt sum( (L-Ln)^2 )
%                       where L: is the length at t_(n+1)
%                             Ln:is the length at t_(n)
%                             dt: time step

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

EnergyBar=0;


%% Loop over Cells
%     % Analytical residual g and Jacobian K
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    nu=Set.nu_Local_EdgeBased/Set.dt;
    for e=1:length(Cell.Ln{i})
        nY=Cell.Cv{i}(e,:);
        Y1=Y.DataRow(nY(1),:);
        if  nY(2) > 0
            Y2=Y.DataRow(nY(2),:);
        else
            Y2=Cell.FaceCentres.DataRow(abs(nY(2)),:);
            nY(2)=abs(nY(2))+Set.NumMainV;
        end
        
        Ln=Cell.Ln{i}(e);
        L=norm(Y1-Y2);
        eij=(Y1'-Y2')./L;
        
        if Set.LocalViscosityOption ==1
            ge=(nu/Ln)* ((L-Ln)/Ln) * [eij ; -eij];
        elseif Set.LocalViscosityOption ==2
            ge=nu* (L-Ln) * [eij ; -eij];
        end
        g=Assembleg(g,ge,nY);
        if nargout>1
            if Set.LocalViscosityOption ==1
                Kij= (nu/Ln)* ((L-Ln)/Ln) * (1/L) * ( eye(3) - eij*eij')...
                    + (nu/Ln^2) * (eij*eij');
            elseif Set.LocalViscosityOption ==2
                Kij= nu* (L-Ln) * (1/L) * ( eye(3) - eij*eij')...
                    + nu * (eij*eij');
            end
            Ke=[ Kij -Kij;
                -Kij  Kij];
            if Set.Sparse
                [si,sj,sv,sk]= KeAssembleSparse(Ke,nY,si,sj,sv,sk);
            else
                K= AssembleK(K,Ke,nY);
            end
            EnergyBar=EnergyBar+ nu/2 *((L-Ln)/Ln)^2;
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end


%%
function K= KeAssemble(K,Ke,nY)

% Assembles Jacobian of a bar of 2 vertices (3x3 components)
dim=3;

for I=1:length(nY) % loop on 3 vertices of triangle
    idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
    idofl=(I-1)*dim+1:I*dim;
    for J=1:length(nY)
        jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
        jdofl=(J-1)*dim+1:J*dim;
        K(idofg,jdofg)=K(idofg,jdofg)+Ke(idofl,jdofl);
    end
end
end


%%
function [si,sj,sv,sk]= KeAssembleSparse(Ke,nY,si,sj,sv,sk)

% Assembles Jacobian of a bar of 2 vertices (6x6 components)
dim=3;

for I=1:length(nY) % loop on 2 vertices of triangle
    idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
    idofl=(I-1)*dim+1:I*dim;
    for J=1:length(nY)
        jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
        jdofl=(J-1)*dim+1:J*dim;
        for d=1:dim
            si(sk+1:sk+dim)=idofg;
            sj(sk+1:sk+dim)=jdofg(d);
            sv(sk+1:sk+dim)=Ke(idofl,jdofl(d));
            sk=sk+dim;
        end
    end
end

end
%%



