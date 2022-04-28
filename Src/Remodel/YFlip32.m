function [Ynew, Tnew] = YFlip32(Ys, Ts, YsToChange, Geo)
	oV = YsToChange;
	n=intersect(intersect(Ts(oV(1),:),Ts(oV(2),:)),Ts(oV(3),:));
	N=unique(Ts(oV,:)); % all nodes
	N=N(~ismember(N,n));

	N3=N(~ismember(N,Ts(oV(1),:)));
	Tnew1=Ts(oV(1),:); Tnew2=Tnew1;
	Tnew1(ismember(Ts(oV(1),:),n(2)))=N3;
	Tnew2(ismember(Ts(oV(1),:),n(1)))=N3;
	Tnew=[Tnew1;Tnew2];

	% The new vertices 
	Xs = zeros(length(n),3);
	for ni = 1:length(n)
		Xs(ni,:) = Geo.Cells(n(ni)).X;
	end
	Ynew=DoFlip32(Ys(oV,:),Xs);
end