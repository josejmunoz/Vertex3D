
function [Si,Sj,Sv]=sReduction(si,sj,sv,sk,Cell)
Sk=0;
for i=1:Cell.n
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end
    Sk=Sk+sk{i};
end 
Si=zeros(Sk,1); % Each vertex affecting 6 nodes
Sj=Si;
Sv=Si;
Sk=0;
for i=1:length(si)
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end
    Si(Sk+1:Sk+sk{i})=si{i}(1:sk{i});
    Sj(Sk+1:Sk+sk{i})=sj{i}(1:sk{i});
    Sv(Sk+1:Sk+sk{i})=sv{i}(1:sk{i});
    Sk=Sk+sk{i};
end 
end 