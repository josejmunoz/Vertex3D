function [Ynew, Tnew] = YFlip32(Ys, Ts, YsToChange, Geo)
	n=intersect(intersect(Ts(YsToChange(1),:),Ts(YsToChange(2),:)),Ts(YsToChange(3),:));
	N=unique(Ts(YsToChange,:)); % all nodes
	N=N(~ismember(N,n));

	N3=N(~ismember(N,Ts(YsToChange(1),:)));
	Tnew1=Ts(YsToChange(1),:); Tnew2=Tnew1;
	Tnew1(ismember(Ts(YsToChange(1),:),n(2)))=N3;
	Tnew2(ismember(Ts(YsToChange(1),:),n(1)))=N3;
	Tnew=[Tnew1;Tnew2];

	% The new vertices 
	Xs = zeros(length(n),3);
	for ni = 1:length(n)
		Xs(ni,:) = Geo.Cells(n(ni)).X;
	end
	Ynew=DoFlip32(Ys(YsToChange,:),Xs);
end