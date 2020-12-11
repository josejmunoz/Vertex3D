function   g=Assembleg(g,ge,nY)
% Assembly of the residual of an element (e.g. Triangle ->length(nY)*3=length(ge)=9,
%                                                edge->    length(nY)*2=length(ge)=6))
dim=3;
for I=1:length(nY) % loop on number of vertices of triangle
    idofg=(nY(I)-1)*dim+1:nY(I)*dim;% global dof
    idofl=(I-1)*dim+1:I*dim;
    g(idofg)=g(idofg)+ge(idofl);
end
end