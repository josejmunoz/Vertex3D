function [g, Energy, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set)
%INITIALIZEKG Summary of this function goes here
%   Detailed explanation goes here
ncell=Cell.n;
% First: Vertices, Faces and Tetrahedrons
dimg=(Set.NumTotalV + Set.NumXs)*3;

if Set.Sparse > 0
    g = sparse(dimg, 1); % Local cell residual
    if nargout > 3
        K = sparse(zeros(dimg));
        if Set.Sparse == 2
            sk=0;
            si=zeros(round((dimg^2)/50),1); % Each vertex is shared by at least 3 cells
            sj=si;
            sv=si;
        end
    end
else
    g = zeros(dimg, 1);
    K = zeros(dimg, dimg);
end

Energy=0;
end

