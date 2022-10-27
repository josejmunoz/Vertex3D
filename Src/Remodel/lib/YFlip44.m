function [Ynew, Tnew] = YFlip44(Ys, Ts, YsToChange, Face, Geo)
    side=[1 2; 2 3; 3 4; 1 4];
	L(1)=norm(Ys(YsToChange(1),:)-Ys(YsToChange(2),:));
	L(2)=norm(Ys(YsToChange(2),:)-Ys(YsToChange(3),:));
	L(3)=norm(Ys(YsToChange(3),:)-Ys(YsToChange(4),:));
	L(4)=norm(Ys(YsToChange(1),:)-Ys(YsToChange(4),:));
	[~,Jun]=min(L);
	VJ=YsToChange(side(Jun,:));
	cVJ3=intersect(Ts(VJ(1),:),Ts(VJ(2),:));
	% cVJ3 is equal
	N=unique(Ts(VJ,:)); % all nodes
	NZ=N(~ismember(N,cVJ3));
	NX=Face.ij;
	Tnew1=Ts(YsToChange(1),:);       Tnew1(ismember(Tnew1,NX(1)))=NZ(~ismember(NZ,Tnew1));
	Tnew2=Ts(YsToChange(2),:);       Tnew2(ismember(Tnew2,NX(2)))=NZ(~ismember(NZ,Tnew2));
	Tnew3=Ts(YsToChange(3),:);       Tnew3(ismember(Tnew3,NX(1)))=NZ(~ismember(NZ,Tnew3));
	Tnew4=Ts(YsToChange(4),:);       Tnew4(ismember(Tnew4,NX(2)))=NZ(~ismember(NZ,Tnew4));
	Tnew=[Tnew1;Tnew2;Tnew3;Tnew4];
	Ynew=DoFlip44(Ys(YsToChange,:),Tnew,L,Geo);
end