function [Geo, g,K,Energy, Set, gr, dyr, dy] = NewtonRaphson(Geo_0, Geo_n, Geo, Dofs, Set, K, g, numStep, t)
	% TODO FIXME There should be a cleaner way for this...
	if Geo.Remodelling
    	dof=Dofs.Remodel;
	else
    	dof=Dofs.Free;
	end
	dy=zeros((Geo.numY+Geo.numF+Geo.nCells)*3, 1);
	dyr=norm(dy(dof)); gr=norm(g(dof));
	gr0=gr;

	fprintf('Step: %i,Iter: %i ||gr||= %e ||dyr||= %e dt/dt0=%.3g\n',numStep,0,gr,dyr,Set.dt/Set.dt0);

	Energy = 0;
    
	Set.iter=1;
    auxgr=zeros(3,1);
    auxgr(1)=gr;
	ig = 1;
	while (gr>Set.tol || dyr>Set.tol) && Set.iter<Set.MaxIter

    	dy(dof)=-K(dof,dof)\g(dof);
    	alpha = LineSearch(Geo_0, Geo_n, Geo, Dofs, Set, g, dy);
    	%% Update mechanical nodes
    	dy_reshaped = reshape(dy * alpha, 3, (Geo.numF+Geo.numY+Geo.nCells))';
    	Geo = UpdateVertices(Geo, Set, dy_reshaped);
		Geo = UpdateMeasures(Geo);

    	%% ----------- Compute K, g ---------------------------------------
    	[g, K, Energy, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
    	dyr=norm(dy(dof)); gr=norm(g(dof));
    	fprintf('Step: % i,Iter: %i, Time: %g ||gr||= %.3e ||dyr||= %.3e alpha= %.3e  nu/nu0=%.3g \n',numStep,Set.iter,t,gr,dyr,alpha,Set.nu/Set.nu0);

    	Set.iter=Set.iter+1;
		auxgr(ig+1)=gr;
		% TODO FIXME, what even is this ?! PVM: In other words, WTF!?
    	if ig ==2
        	ig=0;
    	else
        	ig=ig+1;
    	end
    	if (abs(auxgr(1)-auxgr(2))/auxgr(1)<1e-3 &&...
            	abs(auxgr(1)-auxgr(3))/auxgr(1)<1e-3 &&...
            	abs(auxgr(3)-auxgr(2))/auxgr(3)<1e-3)...
            	|| abs((gr0-gr)./gr0)>1e3
        	Set.iter=Set.MaxIter;
    	end
	end
end

