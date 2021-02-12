function [Cell, Y, Yn, SCn, T, X, Faces, Dofs, Set, Vnew] = backToPreviousStep(Cellp, Yp, Ynp, SCnp, Tp, Xp, Facesp, Dofsp, Setp, Vnewp)
%BACKTOPREVIOUSSTEP Summary of this function goes here
%   Detailed explanation goes here
Cell=Cellp;
Y=Yp;
Yn=Ynp;
SCn=SCnp;
T=Tp;
X=Xp;
Faces=Facesp;
Dofs=Dofsp;
Set=Setp;
Vnew=Vnewp;
Set.N_Rejected_Transfromation=Set.N_Rejected_Transfromation+1;
end

