function [N]=GetN(T,Twb,Xf)


N=(1/4).*ones(size(T));
if nargin > 1 && false 
    ID=1:size(Twb,1);
    IDf=ID(any(ismember(Twb,Xf),2));
    for i=1:length(IDf)
        N(IDf(i),:)=0;
        N(IDf(i),~ismember(T(IDf(i),:),Xf))=[1 1 1]./3;
    end
end 
    



end 