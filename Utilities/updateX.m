function [X]=updateX(X,Y,Ytn,T,Set)


dY=Y.DataOrdered-Ytn(1:Set.NumMainV,:);
nC=zeros(size(X,1),1);
dX=zeros(size(X));
for t=1:size(T.DataOrdered,1)
    if T.NotEmpty(t)
      for n=1:4
         dX(T.DataRow(t,n),:)=dX(T.DataRow(t,n),:)+dY(t,:);
         nC(T.DataRow(t,n))=nC(T.DataRow(t,n))+1;
      end 
    end 
end 

dX(nC>0,:)=dX(nC>0,:)./nC(nC>0);
X=X+dX;