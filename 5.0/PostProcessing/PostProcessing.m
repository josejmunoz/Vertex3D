function PostProcessing(X,Y,Cn,Cv,Ld,Lv,Cell,file,TimeStep)
if TimeStep ==0
    InitVtk(file)
end 
CreateVtkBar(X,Cn,Ld,file,'n',TimeStep)
CreateVtkBar(Y,Cv,Lv,file,'v',TimeStep)
CreateVtkVol(Y,Cell,X,file,TimeStep)
