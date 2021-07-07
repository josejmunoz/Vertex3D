function [g, K, Energy] = KgBulkElem(x, x0, mu, lambda, Neo)
%KGBULKELEM Computes elemental residual and Jacobian for bulk viscoelastic domain
%   Can use:
%   Neo=0: St. Venant-Kirchoff elastic poential:
%   W(E)=lambda tr(E)^2 + 2*mu tr(E^2)
%
%   Neo=1: Neo-Hookean elastic poential:
%   W(E)=lambda ln(J)^2 + mu (tr(E)-ln(J))
%
%   Neo=0: Modififed Neo-Hookean for handling negative volumes (may happen in iterative process)
%   W(E)=lambda ln(J^2)^2 + mu (tr(E)-ln(J^2))
%   If no lambda is given, assumed to be zero (no volumetric stiffness)
%   
%   INPUT:
%   x(i,:)=coordinates of node i (assumed 4 nodes) on deformed element
%   x0(i,:)=coordinates of node i (assumed 4 nodes) on undeformed element
%   mu,lambda = material shear stiffness
% 
%   OUTPUT:
%   g=elemental residual
%   K=elemental Jacobian
%   S=2nf Piola-Kirchhof stress tensor on last GP
%
%   Designed by Jose J. Muñoz

NeoH=2;
if exist('Neo','var')
    NeoH=Neo;
end

if ~exist('lambda','var')
    lambda=0;
end

if size(x,1)~=4
    error('Bulk stiffness can only be used with 4 nodal tetrahedra.')
end

[nnod,dim]=size(x);

ng=1;% Number of Gauss points

if ng==1
    xg=[1 1 1]/4;
    wg=1;
elseif ng==4
    a=0.5854102;
    b=(1-a)/3;
    xg=[ b b b
        a b b
        b a b
        b b a];
    wg=[1 1 1 1]/4;
end

wg=wg/6; % 1/6 because volume of reference tet =1/6

% Derivatives are constant
DN=[1 0 0 -1  % dN/dxi
    0 1 0 -1  % dN/deta
    0 0 1 -1];% dN/dchi

g=zeros(nnod*dim,1);
K=zeros(nnod*dim);
Energy=0;

for ig=1:ng
    % Unnecessary: dxdxi=x'*DN'; % 3x3 Jacobain dx_i/dxi_j
    dXdxi=x0'*DN'; % 3x3     Jacobain dX_i/dxi_j
    gradXN=(dXdxi')\DN; % each column: gradient of each shape function
    F=x'*gradXN';
    gradxN=(F')\gradXN; % each column: gradient of each shape function
    J=det(F);
    E=0.5*(F'*F-eye(3));
    trE=sum(diag(E));
    Je=det(dXdxi);
    lJ=log(J);
    lJ2=log(J^2);
    
    if J<0
%         warning('Inverted tetrahedra');
%         ME = MException('KgBulkElem:invertedTetrahedralElement', ...
%             'Inverted Tetrahedral Element');
%         throw(ME)
    end
    
    if NeoH ~= 2 && Je<0
        error('Tetrahedral Element orientation need to be swapped')
    end
    
    if NeoH==1
        Energy=Energy+0.5*lambda*lJ^2+mu*(trE-lJ)*Je*wg(ig);
    elseif NeoH==2
        Energy=Energy+0.5*lambda*lJ2^2+mu*(trE-lJ2)*Je*wg(ig);
    else % St Venant
        Energy=Energy+0.5*lambda*trE^2+mu*sum(diag(E*E))*Je*wg(ig);
    end
    for a=1:nnod
        gradXNa=gradXN(:,a);
        gradxNa=gradxN(:,a);
        idof=(a-1)*dim+1:a*dim;
        if NeoH==1
            K1=(lambda*lJ-mu)*gradxNa+mu*F*gradXNa;
        elseif NeoH==2
            K1=2*(lambda*lJ2-mu)*gradxNa+mu*F*gradXNa;
        else
            K1=(lambda*trE*F+2*mu*F*E)*gradXNa;
        end
        g(idof)=g(idof)+wg(ig)*K1*Je;
        for b=1:nnod
            gradXNb=gradXN(:,b);
            gradxNb=gradxN(:,b);
            jdof=(b-1)*dim+1:b*dim;
            NxaNxb=gradxNa*gradxNb';
            if NeoH==1
                K1=mu*(gradXNa'*gradXNb)*eye(3);
                K2=lambda*NxaNxb;
                K3=-(lambda*lJ-mu)*NxaNxb';
                Kab=K1+K2+K3;
            elseif NeoH==2
                K1=mu*(gradXNa'*gradXNb)*eye(3);
                K2=4*lambda*NxaNxb;
                K3=-2*(lambda*lJ2-mu)*NxaNxb';
                Kab=K1+K2+K3;
            else
                K1=(trE*lambda*(gradXNa'*gradXNb)+2*mu*(gradXNb'*E*gradXNa))*eye(3);
                K2=lambda*F*gradXNa*(F*gradXNb)'+mu*(F*gradXNb*(F*gradXNa)');
                K3=mu*F*(F')*(gradXNa'*gradXNb);
                Kab=K1+K2+K3;
            end
            K(idof,jdof)=K(idof,jdof)+wg(ig)*Je*Kab;
        end
    end
end
end
