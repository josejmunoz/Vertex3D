function [XgID,X]=SeedWithBoundingBox(X,s)
    % This funcrion seeds nodes in undesired entities (edges, faces and tetrahedrons) 
    % while cell-centers are bounded by ghost nodes. 
    
    nCells = size(X,1); 
    
    r0=mean(X);
    r=5*max(max(abs(X-r0)));

    %% 1)  Define bounding Nodes
    %  Bounding Shpere 
    theta=linspace(0,2*pi,5);
    phi=linspace(0,pi,5);
    [theta,phi]=meshgrid(theta,phi);
    x=r*sin(phi).*cos(theta);
    y=r*sin(phi).*sin(theta);
    z=r*cos(phi);
    x=reshape(x,size(x,1)*size(x,2),1);
    y=reshape(y,size(y,1)*size(y,2),1);
    z=reshape(z,size(z,1)*size(z,2),1);
    Xg=[x y z] + r0;  
    Xg=uniquetol(Xg,'ByRows',1e-6); % FIXME, is this a good idea?
    
    %% 2) Do first Delaunay with ghost nodes
    XgID=size(X,1)+1:size(X,1)+size(Xg,1);
    XgIDBB = XgID;
    X=[X;Xg];
    Twg=delaunay(X);
    
    %% 3) intitilize 
    Side=[1 2 3; 1 2 4; 2 3 4; 1 3 4];
    Edges=[1 2; 2 3; 3 4; 1 3; 1 4; 3 4; 1 4];  
    % find real tests 
    Vol=zeros(size(Twg,1),1);
    AreaFaces=zeros(size(Twg,1)*3,4);
    LengthEdges=zeros(size(Twg,1)*3,6);
    % Volc=0;
    Arc=0;
    Lnc=0;
    
    %%  4) compute the size of Real Entities (edges, faces and tetrahedrons) 
    for i=1:size(Twg,1)  
        %-----------Volume 
    %     if sum(ismember(Twg(i,:),XgID))==0
    %         Vol(i)=abs(    1/6*dot(   cross( X(Twg(i,2),:)-X(Twg(i,1),:),X(Twg(i,3),:)-X(Twg(i,1),:) )  ,...
    %                                                                                      X(Twg(i,4),:)-X(Twg(i,1),:)  )  );
    %         Volc=Volc+1;
    %     else
    %         Vol(i)=0;
    %     end
    
    %    %----------- Area
        for j=1:4
            if sum(ismember(Twg(i,Side(j,:)),XgID))==0
               AreaFaces(i,j)=AreTri(X(Twg(i,Side(j,1)),:),X(Twg(i,Side(j,2)),:),X(Twg(i,Side(j,3)),:));
               Arc=Arc+1;
            else 
               AreaFaces(i,j)=0;
            end 
        end 
        %-----------Length
        for j=1:6
            if sum(ismember(Twg(i,Edges(j,:)),XgID))==0
               LengthEdges(i,j)=norm(X(Twg(i,Edges(j,1)),:)-X(Twg(i,Edges(j,2)),:));
               Lnc=Lnc+1;
            else 
               LengthEdges(i,j)=0;
            end 
        end    
    end 
    
    % mVol=sum(Vol)/Arc;
    % mArea=sum(sum(AreaFaces))/Arc;
    % mLength=sum(sum(LengthEdges))/Lnc;
   
    %% 5) seed nodes in big Entities (based on characteristic Length h) 
    for i=1:size(Twg,1)  
          %---- Seed according to volume  
    %     if sum(ismember(Twg(i,:),XgID))==0 && Vol(i)>VolTol
    %             [X,XgID]=SeedNodeTet(X,XgID,Twg(i,:),h);   
    %     end 
        
        %---- Seed according to area 
        for j=1:4
            if sum(ismember(Twg(i,Side(j,:)),XgID))==0
                if  AreaFaces(i,j)>(s)^2
                    [X,XgID]=SeedNodeTri(X,XgID,Twg(i,Side(j,:)),s);
                end 
            end 
        end  
        
        %---- Seed according to length
        for j=1:6
            if sum(ismember(Twg(i,Edges(j,:)),XgID))==0 && LengthEdges(i,j)>2*s % LengthEdges(i,j)>LengthTol*mLength
    %             [X,XgID]=SeedNodeBar(X,XgID,Twg(i,Edges(j,:)),h);
                [X,XgID]=SeedNodeTet(X,XgID,Twg(i,:),s); 
                break
            end 
        end 
    end 
    
    %% 6)  Seed on ghost Tets
    for i=1:length(Vol)   
        if sum(ismember(Twg(i,:),XgID))>0 
            [X,XgID]=SeedNodeTet(X,XgID,Twg(i,:),s);
         end 
    end

    %% 7)  Remove bounding box nodes
    X(XgIDBB,:) = []; 
    XgID = (nCells+1):size(X,1); 
end


function [X,XgID]=SeedNodeTet(X,XgID,Twgi,h)
    XTet=X(Twgi,:);
    Center=1/4*(sum(XTet,1));
    nX=zeros(4,3);
    for i=1:4
        vc=Center-XTet(i,:);
        dis=norm(vc);
        dir=vc/dis;
        offset=h*dir;
        if dis>norm(offset)
            % offset
            nX(i,:)=XTet(i,:)+offset;
        else 
            % barycenteric
            nX(i,:)=XTet(i,:)+vc;
        end      
    end 
    
    
    nX(ismember(Twgi,XgID),:)=[];
    nX=uniquetol(nX,1e-12*h,'ByRows',true);
    [nX]=CheckReplicateedNodes(X,nX,h);
    nXgID=size(X,1)+1:size(X,1)+size(nX,1);
    X=[X;nX];
    XgID=[XgID nXgID ];
end 


function  [X,XgID]=SeedNodeTri(X,XgID,Tri,h)
    XTri=X(Tri,:);
    Center=1/3*(sum(XTri,1));
    nX=zeros(3,3);
    for i=1:3
        vc=Center-XTri(i,:);
        dis=norm(vc);
        dir=vc/dis;
        offset=h*dir;
        if dis>norm(offset)
            % offset
            nX(i,:)=XTri(i,:)+offset;
        else 
            % barycenteric
            nX(i,:)=XTri(i,:)+vc;
        end      
    end 
    
    nX(ismember(Tri,XgID),:)=[];
    nX=uniquetol(nX,1e-12*h,'ByRows',true);
    [nX]=CheckReplicateedNodes(X,nX,h);
    nXgID=size(X,1)+1:size(X,1)+size(nX,1);
    X=[X;nX];
    XgID=[XgID nXgID ];
end 

function  [X,XgID]=SeedNodeBar(X,XgID,Edge,h)
    XEdge=X(Edge,:);
    Center=1/2*(sum(XEdge,1));
    nX=zeros(2,3);
    for i=1:2
        vc=Center-XEdge(i,:);
        dis=norm(vc);
        dir=vc/dis;
        offset=h*dir;
        if dis>norm(offset)
            % offset
            nX(i,:)=XEdge(i,:)+offset;
        else 
            % barycenteric
            nX(i,:)=XEdge(i,:)+vc;
        end      
    end 
    nX=unique(nX,'row');
    nXgID=size(X,1)+1:size(X,1)+size(nX,1);
    X=[X;nX];
    XgID=[XgID nXgID ]; Main %???
end 

function [nX]=CheckReplicateedNodes(X,nX,h)
    ToBeRemoved=false(size(nX,1),1);
    for jj=1:size(nX,1)
        m=[X(:,1)-nX(jj,1) X(:,2)-nX(jj,2) X(:,3)-nX(jj,3)];
        m=m(:,1).^2+m(:,2).^2+m(:,3).^2;
        m=m.^(1/2);
        m=min(m);
        if m<1e-2*h
            ToBeRemoved(jj)=true;
        end 
    end 
    nX(ToBeRemoved,:)=[];
end 

function [Area]=AreTri(P1,P2,P3)
    Area =1/2*norm(cross(P2-P1,P3-P1));
end 



























