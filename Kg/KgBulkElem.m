function [g, K]=KgBulkElem(x, x0, mu, lambda)
%KGBULKELEM Computes elemental residual and Jacobian for bulk viscoelastic domain
%   Assumes St. Venant-Kirchoff elastic poential:
%   W(E)=lambda tr(E)^2 + 2*mu tr(E^2)
%   If no lmbda is given, assumed to be zero (no volumetric stiffness)
%   
%   INPUT:
%   x(i,:)=coordinates of node i (assumed 4 nodes) on deformed element
%   x0(i,:)=coordinates of node i (assumed 4 nodes) on undeformed element
%   mu,lambda = material shear stiffness
% 
%   OUTPUT:
%   g=eleemntal residual
%   K=elemental Jacobian
%   S=2nf Piola-Kirchhof stress tensor on last GP
%
%   Designed by Jose J. Muñoz

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

for ig=1:ng
    dxdxi=x'*DN'; % 3x3 Jacobain dx_i/dxi_j
    dXdxi=x0'*DN'; % 3x3     Jacobain dX_i/dxi_j
    F=dxdxi/dXdxi;
    E=0.5*(F'*F-eye(3));
    trE=sum(diag(E));
    detJ=det(dXdxi)*det(F); % | dx/dxi
    
    if detJ<0
        error('Tetrahedral Element orientation need to be swapped')
    end
    
    W=0.5*lambda*trE^2+mu*sum(diag(E*E));
    gradXN=(dXdxi')\DN; % each column: gradient of each shape function
    for a=1:nnod
        gradXNa=gradXN(:,a);
        idof=(a-1)*dim+1:a*dim;
        g(idof)=g(idof)+wg(ig)*(W*eye(3)+lambda*trE*F+2*mu*F*E)*gradXNa*detJ;
        for b=1:nnod
            gradXNb=gradXN(:,b);
            jdof=(b-1)*dim+1:b*dim;
            K1=(trE*lambda*(gradXNa'*gradXNb)+2*mu*(gradXNb'*E*gradXNa))*eye(3);
            K2=lambda*F*gradXNa*(F*gradXNb)'+mu*(F*gradXNb*(F*gradXNa)');
            K3=mu*F*(F')*(gradXNa'*gradXNb);
            K4=(lambda*trE*F+2*mu*F*E+W*eye(3))*(gradXNa*gradXNb');
            K5=(gradXNa*gradXNb')*(lambda*trE*F'+2*mu*E*F');
            K(idof,jdof)=K(idof,jdof)+wg(ig)*detJ*(K1+K2+K3+K4+K5);
        end
    end
% Stresses: S=lambda*trE*eye(3)+2*mu*E;
end
end
