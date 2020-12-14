function PostProcessingVTK(X,Y,T,Cn,Cell,file,TimeStep,Set)
% Create VTK files 



CreateVtkVol(Y.DataOrdered,Cell,X,file,TimeStep)



Ln.L=ones(size(Cn,1),1);
Ln.L0=ones(size(Cn,1),1);
CreateVtkBar(X,Cn,Ln,file,'n',TimeStep)
if ~isempty(T)
    CreateVtkTet(X,T,file,TimeStep)
end 

if Set.Confinement, CreateVtkConfinement(Set,file,TimeStep); end 



% if  Set.gVTK 
%     PlotVectorVTK(Y.DataOrdered,Cell,gs,'gsVTK',file,Set.iIncr)
%     PlotVectorVTK(Y.DataOrdered,Cell,gv,'gvVTK',file,Set.iIncr)
%     PlotVectorVTK(Y.DataOrdered,Cell,gf,'gfVTK',file,Set.iIncr)
% end 

end 



% CnNN=Cn(~any(ismember(Cn,XgID),2),:);
% CnNB=Cn(~(sum(ismember(Cn,XgID),2)==2),:);
% CreateVtkBar(X,CnNN,Ln,file,'nNN',TimeStep)
% CreateVtkBar(X,CnNB,Ln,file,'nNB',TimeStep)

% Lv.L=ones(size(Cv,1),1);
% Lv.L0=ones(size(Cv,1),1);
% CreateVtkBar(Y,Cv,Lv,file,'v',TimeStep)