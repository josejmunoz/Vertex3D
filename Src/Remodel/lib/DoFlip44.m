function Yn=DoFlip44(Y,Tnew,L,Geo)
	center=sum(Y,1)./4;
	% L=mean(L)/2;
	L=mean(L);
	
    % TODO FIXME, other way for this???
    c1 = zeros(1,3);
    c2 = zeros(1,3);
    c3 = zeros(1,3);
    c4 = zeros(1,3);
    for t = 1:size(Tnew,2)
        c1=c1 + (Geo.Cells(Tnew(1,t)).X)./4;
        c2=c2 + (Geo.Cells(Tnew(2,t)).X)./4;
        c3=c3 + (Geo.Cells(Tnew(3,t)).X)./4;
        c4=c4 + (Geo.Cells(Tnew(4,t)).X)./4;
    end

% 	c1=sum(Geo.Cells(Tnew(1,:),:).X,1)./4;
% 	c2=sum(Geo.Cells(Tnew(2,:),:).X,1)./4;
% 	c3=sum(Geo.Cells(Tnew(3,:),:).X,1)./4;
% 	c4=sum(Geo.Cells(Tnew(4,:),:).X,1)./4;
	Lc1=c1-center; Lc1=Lc1/norm(Lc1);
	Lc2=c2-center; Lc2=Lc2/norm(Lc2);
	Lc3=c3-center; Lc3=Lc3/norm(Lc3);
	Lc4=c4-center; Lc4=Lc4/norm(Lc4);
	Y1=center+L.*Lc1;
	Y2=center+L.*Lc2;
	Y3=center+L.*Lc3;
	Y4=center+L.*Lc4;
	Yn=[Y1;Y2;Y3;Y4];
end