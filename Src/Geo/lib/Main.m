
% By Jose J. Munoz
function X=Main(dim)
npoints=30;
dim0=2;
if nargin==0
    dim=3;
end
%%
X0=rand(npoints,2);
X2D=X0;
% Rotate to closer (x,y) 2D plane
if dim==3
    Zp=0.0; % perturbation on z
    a=[0.3 0.4 1.0]; %normal of plane containing points =[a(1) a(2) 1], and d distance
    X0=[X0 -a(1)*X0(:,1)-a(2)*X0(:,2)-a(3)-2*Zp*rand(npoints,1)-Zp];
    R=RotationMatrix(X0);
    % plot3(X0(:,1),X0(:,2),X0(:,3),'o')
    X=(R'*X0')';
    X3=X(:,3);
    X2D=X(:,1:2); % Flatten rotated X
    dim0=3;
end
T=delaunay(X2D(:,1),X2D(:,2));
Xf=GetBoundary(T,X2D); % Fixed points
X2D0=X2D;
[X2D,flag,dJ0,dJ]=RegulariseMesh(T,X2D,Xf);
% plot 2D meshes
% initial mesh
Plot2D(dJ,dJ0,T,X2D,X2D0,Xf)
if dim0==3 % Unrotate to original disposition
    X=[X2D X3];
    X=(R*X')';
end
% Plot in 3D
if dim0==3
    Plot3D(dJ,dJ0,T,X,X0);
end
if flag==1
    fprintf('Regularisation successfully finished\n')
else
    fprintf('Regularisation could not be achieved\n')
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
%% 
function R=RotationMatrix(X)
% fit on plane a*x+b*y-z+d=0
x=[X(:,1) X(:,2) ones(size(X,1),1)];
f=zeros(3,1);
A=zeros(3);
for i=1:3
    for j=1:3
        A(i,j)=sum(x(:,i).*x(:,j));
    end
    f(i)=sum(x(:,i).*X(:,3));
end
a=A\f;
% Find rotation of 2D points to be on plane
% plot3(X(:,1),X(:,2),X(:,3),'o');axis equal
n=[-a(1) -a(2) 1]';
n=n/norm(n);
ez=[0 0 1]';
ex=cross(n,ez);
if norm(ex)<1e-6
    ex=[1 0 0]';
end
ex=ex/norm(ex);
thz=acos(ex'*[1 0 0]');
if ex(2)<0
    thz=-thz;
end
thx=acos(n'*ez);
nc=cross(ex,n);
if nc(3)>0
    thx=-thx;
end
Rz=[cos(thz) -sin(thz) 0
    sin(thz) cos(thz) 0
    0          0       1];
Rx=[1 0 0
    0 cos(thx) -sin(thx)
    0 sin(thx) cos(thx)];
R=Rz*Rx;
end
%%
function  Plot2D(dJ,dJ0,T,X2D,X2D0,Xf)
% Plots flat triangulations in 2D
nele=size(T,1);
npoints=size(X2D,1);
figure(1) % initial mesh
clf
for i=1:npoints
    if min(abs(Xf-i))==0
        plot(X2D(i,1),X2D(i,2),'ro')
    else
        plot(X2D(i,1),X2D(i,2),'bo')
    end
    hold on
end
%
for e=1:nele
    Te=[T(e,:) T(e,1)];
    fill(X2D0(Te,1),X2D0(Te,2),dJ0(e))
end
colorbar
title('det(J) Initial')
% Final mesh
figure(2)
clf
plot(X2D(:,1),X2D(:,2),'o')
hold on
for e=1:nele
    Te=[T(e,:) T(e,1)];
    fill(X2D(Te,1),X2D(Te,2),dJ(e))
end
title('det(J) Final')
colorbar
v=version;
if str2double(v(1))<10
    caxis([min(dJ0) max(dJ0)]);
else
    clim([min(dJ0) max(dJ0)]);
end
end
%%
function Plot3D(dJ,dJ0,T,X,X0)
nele=size(T,1);
figure(3)
clf
figure(4)
clf
for e=1:nele
    Te=[T(e,:) T(e,1)];
    figure(3)
    fill3(X0(Te,1),X0(Te,2),X0(Te,3),dJ0(e))
    hold on;
    figure(4)
    fill3(X(Te,1),X(Te,2),X(Te,3),dJ(e))
    hold on;
end
figure(3)
title('det(J) Initial')
axis equal
colorbar
v=version;
if str2double(v(1))<10
    caxis([min(dJ0) max(dJ0)]);
else
    clim([min(dJ0) max(dJ0)]);
end
figure(4)
title('det(J) Final')
axis equal
colorbar
if str2double(v(1))<10
    caxis([min(dJ0) max(dJ0)]);
else
    clim([min(dJ0) max(dJ0)]);
end
end