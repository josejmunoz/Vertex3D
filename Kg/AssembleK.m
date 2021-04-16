function K= AssembleK(K,Ke,nY)
% Assembly of the Jacobian of an element (e.g. Triangle ->length(nY)=3 size(Ke)=9x9,
%                                              edge->   length(nY)=2  size(Ke)=6x6)
dim=3;
K_new = sparse(size(K, 1), size(K, 2));
for I=1:length(nY) % loop on 3 vertices of triangle
    for J=1:length(nY)
        idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
        idofl=(I-1)*dim+1:I*dim;
        jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
        jdofl=(J-1)*dim+1:J*dim;
        K_new(idofg,jdofg) = Ke(idofl,jdofl);
    end
end
K = K + K_new;
end
