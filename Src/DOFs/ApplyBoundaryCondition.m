function [Geo, Dofs] = applyBoundaryCondition(t, Geo, Dofs, Set)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ApplyBoundaryCondition:											
	%   Modify the DOFs object by including the prescribed points for a
	%   given time. Also update positions of the points in Geo object. 
	% 	Dofs contains the UNVECTORIZED ids (every point has            
	%   a single id only) for the constrained FixC (no displacement    
	%   throughout the simulation) and the prescribed FixP (values may 
	%	displacement may vary throughout the simulation					
	% Input:															
	%   t    : Current time of simulation								
	%   Geo  : Geometry object										   
	%   Dofs : Dofs object. Contains Free, FixP, FixC and Fix			
	%   Set  : User input settings struct                              
	% Output:															
	%   Geo  : Geometry object with updated positions                  
	%   Dofs : DOFS object with updated degrees of freedom				
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if t<=Set.TStopBC && t>=Set.TStartBC
		[dimP, FixIDs] = ind2sub([3, Geo.numY+Geo.numF+Geo.nCells],Dofs.FixP);
		if Set.BC==1
			Geo = UpdateDOFsStretch(FixIDs, Geo, Set);
		elseif Set.BC==2
			[Geo, Dofs] = UpdateDOFsCompress(Geo, Set);
		end
		Dofs.Free(ismember(Dofs.Free,Dofs.FixP))=[];
		Dofs.Free(ismember(Dofs.Free,Dofs.FixC))=[];
	elseif Set.BC==1 || Set.BC==2
		%Dofs.Free=unique([Dofs.Free]);
    end
    
		
end

