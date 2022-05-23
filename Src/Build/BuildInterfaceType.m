function ftype = BuildInterfaceType(ij, XgID, XgTop, XgBottom)
    valueset = 0:2;
    catnames = {'Top' 'Cell-Cell' 'Bottom'};
    categoricalValues = categorical(0:2, valueset,catnames,'Ordinal',true);

	if any(ismember(ij, XgID)) %External: either top or bottom
        if any(ismember(ij, XgTop)) %Top
            ftype = categoricalValues(1);
        elseif any(ismember(ij, XgBottom)) % Bottom/Substrate'Top'
            ftype = categoricalValues(3);
        else % Border face
            ftype = categoricalValues(2);
        end
    else % Lateral domain/cell-cell contact
		ftype = categoricalValues(2);
	end
end