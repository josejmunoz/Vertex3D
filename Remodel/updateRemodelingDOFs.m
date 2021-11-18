function [Dofs] = updateRemodelingDOFs(Dofs, nV, nC, Y)
%UPDATEREMODELINGDOFS Summary of this function goes here
%   Detailed explanation goes here
    if ~isempty(nV)
        nV=reshape(nV,1,length(nV));
        %Create remodeling DOFs
        YDofs=3.*(kron(nV,[1 1 1])-1)+kron(ones(1,length(nV)),[1 2 3]); 
        if ~isempty(nC)
            faceDofs=3.*(kron(nC,[1 1 1])-1)+kron(ones(1,length(nC)),[1 2 3]);
        else 
            faceDofs=[];
        end 
        Dofs.Remodel=[YDofs faceDofs+Y.n*3];
        
        %Dofs.Remodel = setdiff(Dofs.Remodel, Dofs.Debris);
    else
        Dofs.Remodel = [];
    end
end

