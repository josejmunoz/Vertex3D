function g = Assembleg(g,ge,nY)
% Assembly of the residual of an element (e.g. Triangle ->length(nY)*3=length(ge)=9,
%                                                edge->    length(nY)*2=length(ge)=6))
dim=3;
idofg = zeros(length(nY), 1);

for I=1:length(nY) % loop on number of vertices of triangle
    idofg((I-1)*dim+1:I*dim) = (nY(I)-1)*dim+1:nY(I)*dim;% global dof
end

g(idofg)=g(idofg)+ge;
end