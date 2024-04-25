function [Geo, g,K,Energy, Set, gr, dyr, dy] = NewtonRaphson(Geo_0, Geo_n, Geo, Dofs, Set, K, g, numStep, t)
	% TODO FIXME There should be a cleaner way for this...
	if Geo.Remodelling
    	dof=Dofs.Remodel;
	else
    	dof=Dofs.Free;
	end
	dy=zeros((Geo.numY+Geo.numF+Geo.nCells)*3, 1);
	dyr=norm(dy(dof)); 
    gr=norm(g(dof));
	gr0=gr;

	Geo.log = sprintf('%s Step: %i,Iter: %i ||gr||= %e ||dyr||= %e dt/dt0=%.3g\n',Geo.log, numStep,0,gr,dyr,Set.dt/Set.dt0);

	Energy = 0;
    
	Set.iter=1;
    auxgr=zeros(3,1);
    auxgr(1)=gr;
	ig = 1;

	dy(dof)=-Set.dt/Set.nu * g(dof);
    
	%% Update mechanical nodes
	dy_reshaped = reshape(dy, 3, (Geo.numF+Geo.numY+Geo.nCells))';
	Geo = UpdateVertices(Geo, Set, dy_reshaped);
	Geo = UpdateMeasures(Geo);
end

