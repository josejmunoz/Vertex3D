function [X]=Example(e)


if e==1
    X=0;
    Y=0:2;
    [X,Y]=meshgrid(X,Y);
    X=reshape(X,size(X,1)*size(X,2),1);
    Y=reshape(Y,size(Y,1)*size(Y,2),1);
    X=[X Y zeros(length(X),1)+rand(length(X),1)*0];
    
elseif e==2
    %     % Sphere
        [X,Y,Z,~] = mySphere(40);
        X=[X' Y' Z'].*2;
        X(:,3)=X(:,3);
    %     % flat 
        r=2;
        n=6;
        x=linspace(-r-r/2,r+r/2,n);   
        y=linspace(-r-r/2,r+r/2,n);   
        [x,y]=meshgrid(x,y);
        z=-3-.1.*ones(size(x));
        x=reshape(x,size(x,1)*size(x,2),1);
        y=reshape(y,size(y,1)*size(y,2),1);
        z=reshape(z,size(z,1)*size(z,2),1);
        XX=[x y z];
        XX(:,[1 2])=XX(:,[1 2])+0*rand(size(XX(:,[1 2])))./190000;
        X=[X;XX];
elseif e==3
     n=20;
    t = linspace(0,2*pi,n);

    r = ones(n,1)*3;
    [X,Y,Z] = cylinder(r,10);
    X=reshape(X,size(X,1)*size(X,2),1);
    Y=reshape(Y,size(Y,1)*size(Y,2),1);
    Z=reshape(Z,size(Z,1)*size(Z,2),1);

    X=[X Y Z.*n*2];
    X=X(:,:);
    X=unique(X,'rows');    
    k=0;
    for i=1:10
    X(k+1:k+n,1)=X(k+1:k+n,1)+(cos(t)+1)'*6;
    k=k+n;
    end 
    
elseif e==4
    X=0:2;
    Y=0:2;
    [X,Y]=meshgrid(X,Y);
    X=reshape(X,size(X,1)*size(X,2),1);
    Y=reshape(Y,size(Y,1)*size(Y,2),1);
    X=[X Y ones(length(X),1)+rand(length(X),1)*0];
    
elseif e==5
        X=0:2;
    Y=0:2;
    [X,Y]=meshgrid(X,Y);
    X=reshape(X,size(X,1)*size(X,2),1);
    Y=reshape(Y,size(Y,1)*size(Y,2),1);
    X=[X Y ones(length(X),1)+rand(length(X),1)*0;
       X Y ones(length(X),1)*2];
end 