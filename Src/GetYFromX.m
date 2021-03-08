function [Y] = GetYFromX(X,XgID,T,f)
%GETYFROMX Computes vertex positions (Y) form nodal position X
%   Detailed explanation goes heres
dim=size(X,2);
nvert=size(T,1);
Y=zeros(nvert,dim);
for i=1:nvert
    x=X(T(i,:),:); % 3 nodes taht interpolate vertex i
    if abs(sum(ismember(T(i,:),XgID))-3)<eps
        x=X(T(i,:),:);
        Center=1/4*(sum(x,1));
        vc=Center-X(T(i,~ismember(T(i,:),XgID)),:);
        dis=norm(vc);
        dir=vc/dis;
        offset=f*dir;
        Y(i,:)=X(T(i,~ismember(T(i,:),XgID)),:)+offset;
    else 
        for n=1:size(T,2)
             Y(i,1:3)=Y(i,1:3)+(1/4)*x(n,1:3);
        end
    end 
end
end

