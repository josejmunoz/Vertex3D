function [Ynew, Tnew, TRemoved, treeOfPossibilities, arrayPos] = YFlipNM_recursive(TOld, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect, treeOfPossibilities, parentNode, arrayPos)

endNode = 2;

Told_original = TOld;
if size(TOld, 1) == 3
    [Ynew_c, Tnew_c] = YFlip32(oldYs, TOld, [1 2 3], Geo);
    TRemoved{arrayPos} = TOld;
    Tnew{arrayPos} = Tnew_c;
    Ynew{arrayPos} = Ynew_c;
    treeOfPossibilities = addedge(treeOfPossibilities, parentNode, arrayPos);
    treeOfPossibilities = addedge(treeOfPossibilities, arrayPos, endNode);
    arrayPos = arrayPos + 1;
else
    % https://link.springer.com/article/10.1007/s00366-016-0480-z#Fig2
    for numPair = 1:size(possibleEdges, 1)
        [valence, sharedTets, tetIds] = edgeValenceT(Told_original, possibleEdges(numPair, :));
        
        % Valence == 1, is an edge that can be removed.
        % Valence == 2, a face can be removed.
        if valence == 2
            [Ynew_23, Tnew_23] = YFlip23(oldYs, Told_original, tetIds, Geo);

            TRemoved{arrayPos} = Told_original(tetIds, :);
            Tnew{arrayPos} = Tnew_23;
            Ynew{arrayPos} = Ynew_23;
            treeOfPossibilities = addedge(treeOfPossibilities, parentNode, arrayPos);
            

            TOld = Told_original;
            TOld(tetIds, :) = [];
            TOld = vertcat(TOld, Tnew_23);

            %Update and get the tets that are associated to that
            %edgeToDisconnect
            % Valence should have decreased
            [~, TOld_new, ~] = edgeValenceT(TOld, XsToDisconnect);
            [Ynew, Tnew, TRemoved, treeOfPossibilities, arrayPos] = YFlipNM_recursive(TOld_new, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect, treeOfPossibilities, arrayPos, arrayPos+1);
            % What should I do after I receive the Tnew from other
            % How can I get all the correct new tets?
        end
    end
end
end