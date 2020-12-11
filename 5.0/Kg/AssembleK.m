function K= AssembleK(K,Ke,nY)
% Assembly of the Jacobian of an element (e.g. Triangle ->length(nY)=3 size(Ke)=9x9,
%                                              edge->   length(nY)=2  size(Ke)=6x6)
dim=3;
for I=1:length(nY) % loop on 3 vertices of triangle
    for J=1:length(nY)
        idofg=(nY(I)-1)*dim+1:nY(I)*dim; % global dof
        idofl=(I-1)*dim+1:I*dim;
        jdofg=(nY(J)-1)*dim+1:nY(J)*dim; % global dof
        jdofl=(J-1)*dim+1:J*dim;
        K(idofg,jdofg)=K(idofg,jdofg)+Ke(idofl,jdofl);
    end
end
end
