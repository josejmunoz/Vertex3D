function  [X,flag,dJ0,dJ,Xf]=RegulariseMesh(T,X,Xf)
% Moves nodes in order to get regular triangles in a 2D triangular mesh
% T(e,:)=connectivity of triangle e (3 compoenents)
% X(i,:) coordinate of node i
% Xf = list of nodes that will not move. If not given, assumes those nodes at the
% boundary.
% By Jose J. Munoz
if nargin==2
    Xf=GetBoundary(T,X);
end
[x,dofc,doff]=Getx(X,Xf);
dJ0=ComputeJ(T,X);
fun=@(x) gBuild(x,dofc,doff,T,X);
[x,~,flag]=fsolve(fun,x);
[np,dim]=size(X);
xf=reshape(X',[dim*np,1]);
xf(doff)=x;
X=reshape(xf,[dim,np])';
dJ=ComputeJ(T,X);
end
%%
function g=gBuild(x,dofc,doff,T,X)
[nele,npe]=size(T);
[np,dim]=size(X);
xf=reshape(X',[dim*np,1]);
xf(doff)=x;
X=reshape(xf,[dim,np])';
%
xi=sqrt(3)/3*[1 1];
[~,dN]=ShapeFunctions(xi);
g=zeros(np*dim,1);
for e=1:nele
    Xe=X(T(e,:),:);
    Je=Xe'*dN;
    dJ=det(Je);
    fact=(dJ-1)*dJ;
    for i=1:npe
        idof=(T(e,i)-1)*dim+1:T(e,i)*dim;
        g(idof)=g(idof)+fact*(Je')\(dN(i,:)');
    end
end
g(dofc)=[];
end
%
function dJ=ComputeJ(T,X)
nele=size(T,1);
dJ=zeros(nele,1);
xi=sqrt(3)/3*[1 1];
[~,dN]=ShapeFunctions(xi);
for e=1:nele
    Xe=X(T(e,:),:);
    Je=Xe'*dN;
    dJ(e)=det(Je);
end
end
function [x,dofc,doff]=Getx(X,Xf)
[np,dim]=size(X);
x=reshape(X',[dim*np,1]);
doff=1:np*dim;
if nargin>1
    nf=length(Xf);
    dofc=kron((Xf-1)*dim,[1,1])+kron(ones(1,nf),[1,2]);
    doff(dofc)=[];
    x(dofc)=[];
end
end
function [N,dN]=ShapeFunctions(xi)
n=length(xi);
N=[];
if n==2
    N=[1-xi(1)-xi(2) xi(1) xi(2)];
    dN=[-1 -1
         1 0
         0 1];
end
end
%%
function nodesExt=GetBoundary(T,X)
np=size(X,1);
nele=size(T,1);
nodesExt=zeros(1,np);
for e=1:nele
    Te=[T(e,:) T(e,1)];
    Sides=[0 0 0];
    for s=1:3
        n=Te(s:s+1);
        for d=1:nele
            if sum(ismember(n,T(d,:)))==2 && d~=e
                Sides(s)=1;
                break;
            end
        end
        if Sides(s)==0
            nodesExt(Te(s:s+1))=Te(s:s+1);
        end
    end
end
nodesExt(nodesExt==0)=[];
end