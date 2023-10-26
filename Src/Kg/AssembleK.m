function K = AssembleK(K,Ke,nY)
% Assembly of the Jacobian of an element (e.g. Triangle ->length(nY)=3 size(Ke)=9x9,
%                                              edge->   length(nY)=2  size(Ke)=6x6)
dim=3;
idofg = zeros(length(nY)*dim, 1);
jdofg = idofg;

for I = 1:length(nY)
    idofg((I-1)*dim+1:I*dim) = (nY(I)-1)*dim+1:nY(I)*dim;
    jdofg((I-1)*dim+1:I*dim) = (nY(I)-1)*dim+1:nY(I)*dim; % global dof
end

% % Compute global DOF indices in a vectorized manner
% [numI, numJ] = meshgrid(idofg, jdofg);
% idx = sub2ind(size(K), numI(:), numJ(:));
% 
% % Reshape Ke to match the linear indices for assembly
% Ke_flat = reshape(Ke, dim^4, 1);

% Update the matrix K using sparse matrix addition
%K(idx) = K(idx) + Ke_flat;
K(idofg,jdofg) = K(idofg,jdofg) + Ke;

end
