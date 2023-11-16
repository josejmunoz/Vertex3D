function ftype = BuildInterfaceType(ij, XgID, XgTop, XgBottom)
    if any(ismember(ij, XgID)) %External: either top or bottom
        if any(ismember(ij, XgTop)) %Top
            ftype = 1;
        elseif any(ismember(ij, XgBottom)) % Bottom/Substrate
            ftype = 3;
        else % Border face
            ftype = 2;
        end
    else % Lateral domain/cell-cell contact
		ftype = 2;
	end
end