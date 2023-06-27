function [alpha]=LineSearch(Geo_0, Geo_n, Geo, Dofs, Set, gc, dy)
	
	%% Update mechanical nodes
	dy_reshaped = reshape(dy, 3, (Geo.numF+Geo.numY+Geo.nCells))';
	
	[Geo] = UpdateVertices(Geo, Set, dy_reshaped);
	Geo   = UpdateMeasures(Geo);

	g = KgGlobal(Geo_0, Geo_n, Geo, Set);
	dof = Dofs.Free;
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