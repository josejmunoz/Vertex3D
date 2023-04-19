function [Ynew, Tnew, TRemoved] = YFlipNM_recursive(TOld, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect)
if size(TOld, 1) == 3
    [Ynew, Tnew] = YFlip32(oldYs, TOld, [1 2 3], Geo);
    TRemoved = TOld;
else
    % https://link.springer.com/article/10.1007/s00366-016-0480-z#Fig2
    possibleEdgesToKeep = [];
    for numPair = 1:size(possibleEdges, 1)
        [valence, sharedTets, tetIds] = edgeValenceT(TOld, possibleEdges(numPair, :));
        
        % Valence == 1, is an edge that can be removed.
        % Valence == 2, a face can be removed.
        if valence == 2
            possibleEdgesToKeep(end+1, :) = possibleEdges(numPair, :);
            [Ynew, Tnew_23] = YFlip23(oldYs, TOld, tetIds, Geo);
            TOld(tetIds, :) = [];
            TOld = vertcat(TOld, Tnew_23);
            %Update and get the tets that are associated to that
            %edgeToDisconnect
            % Valence should have decreased
            [~, TOld, ~] = edgeValenceT(TOld, XsToDisconnect);
            [Ynew, Tnew, TRemoved] = YFlipNM_recursive(TOld, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect);
            TOld = vertcat(TOld, Tnew);
            % What should I do after I receive the Tnew from other
            % How can I get all the correct new tets?
        end
    end
end
end