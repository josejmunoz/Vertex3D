function [IsNotConvex,itet]=CheckConvexityCondition(Tnew,T)
%CHECKCONVEXITYCONDITION Summary of this function goes here
%   Detailed explanation goes here

IsNotConvex=false;
itet=[];
for i=1:T.n
    Tlogical=ismember(Tnew,T.DataRow(i,:));
    Tlogical=all(Tlogical,2);
    if any(Tlogical)
       IsNotConvex=true;
       itet=i;
       return
    end 
end 

end

