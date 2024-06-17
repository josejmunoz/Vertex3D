function ftype = BuildInterfaceType(ij, XgID, XgTop, XgBottom)
    if any(ismember(ij, XgTop)) %Top
        ftype = 1;
    elseif any(ismember(ij, XgBottom)) % Bottom/Substrate
        ftype = 3;
    else% Border face
        ftype = 2;
    end
end