function Y=GetYSubstrate(Y,X,T,XgSub,XgID,f,S)

nvert=size(Y,1);
for i=1:nvert
    aux=ismember(T(i,:),XgSub);
    if abs(sum(aux))>eps
        XX=X(T(i,~ismember(T(i,:),XgID)),:);
        if size(XX,1)==1
            x=X(T(i,~aux),:);
            Center=1/3*(sum(x,1));
            vc=Center-X(T(i,~ismember(T(i,:),XgID)),:);
            dis=norm(vc);
            dir=vc/dis;
            offset=f*dir;
            Y(i,:)=X(T(i,~ismember(T(i,:),XgID)),:)+offset;
                    Y(i,3)=S;

        elseif size(XX,1)==2
            X12=XX(1,:)-XX(2,:);
            ff=sqrt(f^2-(norm(X12)/2)^2);
            XX=sum(XX,1)/2;
            Center=1/3*(sum(X(T(i,~ismember(T(i,:),XgSub)),:),1));
            vc=Center-XX;
            dis=norm(vc);
            dir=vc/dis;
            offset=ff*dir;
            Y(i,:)=XX+offset;
                    Y(i,3)=S;
        elseif size(XX,1)==3
            Y(i,:)=(1/3).*(sum(X(T(i,~ismember(T(i,:),XgSub)),:),1));
                    Y(i,3)=S;

        end 
    end 
end
end 