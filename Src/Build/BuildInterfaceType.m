function ftype = BuildInterfaceType(ij, XgID, XgTop, XgBottom)
	if any(ismember(ij, XgID)) %External: either top or bottom
        if any(ismember(ij, XgTop)) %Top
            ftype = 0;
        elseif any(ismember(ij, XgBottom)) % Bottom/Substrate
            ftype = 2;
        else % Border face
            ftype = 1;
        end
    else % Lateral domain/cell-cell contact
		ftype = 1;
	end
end