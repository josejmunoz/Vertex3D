function Centre = BuildFaceCentre(ij, ncells, X, Ys, H, extrapolateFaceCentre) 
	% TODO FIXME This function does much more in its original version. For 
	% now, let's only calculate the interpolation. 
	% To add: 
	% - Check if face center (from opposite node) has already been built 
	% - Orientation logic ? 
    Centre=sum(Ys,1)/length(Ys); 
    if sum(ismember(ij,1:ncells))==1 && extrapolateFaceCentre
        runit=(Centre-X);
        runit=runit/norm(runit);
        Centre=X+H.*runit;
    end
end 
 
