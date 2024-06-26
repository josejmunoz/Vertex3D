function [contractilityValue, Geo] = getContractilityBasedOnLocation(currentFace, currentTri, Geo, Set)
%GETCONTRACTILITYBASEDONLOCATION Summary of this function goes here
%   Detailed explanation goes here
    
    CUTOFF = 100;
    
    if isempty(currentTri.ContractilityValue)
        
        if Set.DelayedAdditionalContractility == 1
            contractilityValue = getDelayedContractility(Set.currentT, Set.purseStringStrength, currentTri, CUTOFF * Set.purseStringStrength);
        else
            contractilityValue = getIntensityBasedContractility(Set, currentFace);
        end

        switch (currentFace.InterfaceType)
            case 1 % Top
                if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                    contractilityValue = contractilityValue * Set.cLineTension;
                else
                    contractilityValue = Set.cLineTension;
                end
            case 2 % Lateral
                %% DO LATERAL CABLES HAVE A DIFFERENT MINUTES DELAY
                %% CAN IT BE BASED ON HOW FAST IT IS STRAINED?
                if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                    contractilityValue = contractilityValue * Set.cLineTension;
                else
                    contractilityValue = Set.cLineTension/100;
                end
            case 3 % Bottom/Substrate
                contractilityValue = Set.cLineTension/100;
            otherwise
                contractilityValue = Set.cLineTension;
        end
        %% Adding noise to contractility
        contractilityValue = AddNoiseToParameter(contractilityValue, Set.noiseContractility, currentTri);
    
        for cellToCheck = currentTri.SharedByCells
            facesToCheck = Geo.Cells(cellToCheck).Faces;
            faceToCheckID = ismember(sort(vertcat(facesToCheck.ij), 2), sort(currentFace.ij, 2), "rows");
            if any(faceToCheckID)
                trisToCheck = Geo.Cells(cellToCheck).Faces(faceToCheckID).Tris;
                for triToCheckID = 1:length(trisToCheck)
                    triToCheck = trisToCheck(triToCheckID);
                    if all(ismember(sort(triToCheck.SharedByCells), sort(currentTri.SharedByCells)))
                        Geo.Cells(cellToCheck).Faces(faceToCheckID).Tris(triToCheckID).ContractilityValue = contractilityValue;
                    end
                end
            end
        end
    else
        contractilityValue = currentTri.ContractilityValue;
    end

end

