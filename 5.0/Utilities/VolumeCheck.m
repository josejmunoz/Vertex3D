function IsConsistent=VolumeCheck(Cell,Y)

Sides=[1 2;
       2 3;
       3 1];
IsConsistent=true;
for c=1:Cell.n 
    Tris=Cell.Tris{c};
    Tris(Tris(:,3)>0,3)=Tris(Tris(:,3)>0,3)+Y.n;
    Tris(Tris(:,3)<0,3)=abs(Tris(Tris(:,3)<0,3));
    Id=1:size(Tris,1);
    for t=1:size(Tris,1)
        Side1=Tris(t,Sides(1,:));
        TID=Id(sum(ismember(Tris,Side1),2)==2); TID(TID==t)=[];
        if (Side1(1)==Tris(TID,1) && Side1(2)==Tris(TID,2)) || ...
           (Side1(1)==Tris(TID,2) && Side1(2)==Tris(TID,3)) || ...
           (Side1(1)==Tris(TID,3) && Side1(2)==Tris(TID,1))
           fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
           IsConsistent=false;
           return
        end 
        Side2=Tris(t,Sides(2,:));
        TID=Id(sum(ismember(Tris,Side2),2)==2); TID(TID==t)=[];
        if (Side2(1)==Tris(TID,1) && Side2(2)==Tris(TID,2)) || ...
           (Side2(1)==Tris(TID,2) && Side2(2)==Tris(TID,3)) || ...
           (Side2(1)==Tris(TID,3) && Side2(2)==Tris(TID,1))
           fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
           IsConsistent=false;
           return
        end
        Side3=Tris(t,Sides(3,:));
        TID=Id(sum(ismember(Tris,Side3),2)==2); TID(TID==t)=[];
        if (Side3(1)==Tris(TID,1) && Side3(2)==Tris(TID,2)) || ...
           (Side3(1)==Tris(TID,2) && Side3(2)==Tris(TID,3)) || ...
           (Side3(1)==Tris(TID,3) && Side3(2)==Tris(TID,1))
           fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
           IsConsistent=false;
           return
        end
    end 
end 



end 