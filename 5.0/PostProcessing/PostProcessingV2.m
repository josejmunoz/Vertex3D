function PostProcessingV2(X,Y,Cn,Cv,Cell,file,TimeStep,XgID)
if TimeStep ==0
%     InitVtk(file)
end 

% CnNN=Cn(~any(ismember(Cn,XgID),2),:);
% CnNB=Cn(~(sum(ismember(Cn,XgID),2)==2),:);
% CnBB=Cn(any(ismember(Cn,XgID),2),:);


CreateVtkVol(Y,Cell,X,file,TimeStep)
% Lv.L=ones(size(Cv,1),1);
% Lv.L0=ones(size(Cv,1),1);
Ln.L=ones(size(Cn,1),1);
Ln.L0=ones(size(Cn,1),1);
% CreateVtkBar(Y,Cv,Lv,file,'v',TimeStep)
CreateVtkBar(X,Cn,Ln,file,'n',TimeStep)
% CreateVtkBar(X,CnNN,Ln,file,'nNN',TimeStep)
% CreateVtkBar(X,CnNB,Ln,file,'nNB',TimeStep)