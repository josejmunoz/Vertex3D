% TODO FIXME; this function seems to complex for what it does ?
function Y=BuildYSubstrate(Cell, Cells, XgID, Set, XgSub)
	Tets = Cell.T;
	Y = Cell.Y;
	nverts = length(Tets);
	X = zeros(length(Cells),3);
	for c = 1:length(Cells)
		X(c,:) = Cells(c).X;
	end
	for i=1:nverts
    	aux=ismember(Tets(i,:),XgSub);
    	if abs(sum(aux))>eps
        	XX=X(Tets(i,~ismember(Tets(i,:),XgID)),:);
        	if size(XX,1)==1
            	x=X(Tets(i,~aux),:);
            	Center=1/3*(sum(x,1));
            	vc=Center-X(Tets(i,~ismember(Tets(i,:),XgID)),:);
            	dis=norm(vc);
            	dir=vc/dis;
            	offset=Set.f*dir;
            	Y(i,:)=X(Tets(i,~ismember(Tets(i,:),XgID)),:)+offset;
            	Y(i,3)=Set.SubstrateZ;
        	elseif size(XX,1)==2
            	X12=XX(1,:)-XX(2,:);
            	ff=sqrt(Set.f^2-(norm(X12)/2)^2);
            	XX=sum(XX,1)/2;
            	Center=1/3*(sum(X(Tets(i,~ismember(Tets(i,:),XgSub)),:),1));
            	vc=Center-XX;
            	dis=norm(vc);
            	dir=vc/dis;
            	offset=ff*dir;
            	Y(i,:)=XX+offset;
            	Y(i,3)=Set.SubstrateZ;
        	elseif size(XX,1)==3
            	Y(i,:)=(1/3).*(sum(X(Tets(i,~ismember(Tets(i,:),XgSub)),:),1));
            	Y(i,3)=Set.SubstrateZ;
        	end 
    	end 
	end
end 