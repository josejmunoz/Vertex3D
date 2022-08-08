function [Ynew, Tnew] = YFlip23(Ys, Ts, YsToChange, Geo)
	n3=Ts(YsToChange(1),  ismember(Ts(YsToChange(1),:), Ts(YsToChange(2),:)));
	n1=Ts(YsToChange(1), ~ismember(Ts(YsToChange(1),:),n3));
	n2=Ts(YsToChange(2), ~ismember(Ts(YsToChange(2),:),n3));
	num=[1 2 3 4];
	num=num(Ts(YsToChange(1),:)==n1);
	if num == 2 || num == 4
		Tnew=[n3([1 2]) n2 n1;
			  n3([2 3]) n2 n1;
			  n3([3 1]) n2 n1];
	else
		Tnew=[n3([1 2]) n1 n2;
			  n3([2 3]) n1 n2;
			  n3([3 1]) n1 n2];       
	end
	
	ghostNodes = ismember(Tnew,Geo.XgID);
	ghostNodes = all(ghostNodes,2);
	
	Ynew=DoFlip23(Ys(YsToChange,:),Geo,n3);
	Ynew(ghostNodes,:)=[];
end