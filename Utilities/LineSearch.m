function [alpha]=LineSearch(Cell,Faces, SCn, y0, y,yn,dy,gc,dof,Set,Y,Yn, CellInput)

y0=y;

alpha=1;
y=y0 + alpha*dy; % update nodes

Yt=reshape(y,3,Set.NumTotalV)';
Y=Y.Modify(Yt(1:Set.NumMainV,:));
Cell.FaceCentres=Cell.FaceCentres.Modify(Yt(Set.NumMainV+1:Set.NumTotalV,:));



[g]=KgGlobal(Cell,Faces,SCn, Y0, Y,Yn,y,yn,Set,CellInput);
gr0=norm(gc(dof));   
gr=norm(g(dof));   

if gr0<gr
    R0=dy(dof)'*gc(dof);
    R1=dy(dof)'*g(dof);
    
    R=(R0/R1);
    alpha1=(R/2)+sqrt((R/2)^2-R);
    alpha2=(R/2)-sqrt((R/2)^2-R);
    
    
    if isreal(alpha1) && alpha1<2 && alpha1>1e-3
        alpha=alpha1;
    elseif isreal(alpha2) && alpha2<2 && alpha2>1e-3
        alpha=alpha2;
    else
        alpha=0.1;
    end
else
    alpha=1;
end

end