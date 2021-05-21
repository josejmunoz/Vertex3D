function [g,K]=KgBulk(X,X0,T,mu,lambda)

nele=size(T,1);
[nnod,dim]=size(X);
dof=nnod*dim;
nnode=size(T,2);
g=zeros(dof,1);
K=zeros(dof);
for e=1:nele
    Te=T(e,:);
    edof=kron(ones(nnode,1),(1:dim)')+kron((Te-1)'*dim,ones(dim,1));
    Xe=X(Te,:);
    Xe0=X0(Te,:);
    [ge,Ke]=gKBulkElem(Xe,Xe0,mu,lambda);
    g(edof)=g(edof)+ge;
    K(edof,edof)=K(edof,edof)+Ke;
end
end