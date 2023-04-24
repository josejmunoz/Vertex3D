function [Ynew, Tnew, TRemoved] = YFlipNM_recursive(TOld, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect, treeDepth)
Told_original = TOld;
if size(TOld, 1) == 3
    [Ynew_c, Tnew_c] = YFlip32(oldYs, TOld, [1 2 3], Geo);
    TRemoved{treeDepth} = vertcat(TRemoved{treeDepth}, TOld);
    Tnew{treeDepth} = vertcat(Tnew{treeDepth}, Tnew_c);
    Ynew{treeDepth} = vertcat(Ynew{treeDepth}, Ynew_c);
else
    % https://link.springer.com/article/10.1007/s00366-016-0480-z#Fig2
    possibleEdgesToKeep = [];
    for numPair = 1:size(possibleEdges, 1)
        [valence, sharedTets, tetIds] = edgeValenceT(Told_original, possibleEdges(numPair, :));
        
        % Valence == 1, is an edge that can be removed.
        % Valence == 2, a face can be removed.
        if valence == 2
            possibleEdgesToKeep(end+1, :) = possibleEdges(numPair, :);
            [Ynew_23, Tnew_23] = YFlip23(oldYs, Told_original, tetIds, Geo);

            TRemoved{treeDepth} = vertcat(TRemoved{treeDepth}, Told_original(tetIds, :));
            Tnew{treeDepth} = vertcat(Tnew{treeDepth}, Tnew_23);
            Ynew{treeDepth} = vertcat(Ynew{treeDepth}, Ynew_23);

            TOld = Told_original;
            TOld(tetIds, :) = [];
            TOld = vertcat(TOld, Tnew_23);

            %Update and get the tets that are associated to that
            %edgeToDisconnect
            % Valence should have decreased
            TRemoved{treeDepth+1} = [];
            Tnew{treeDepth+1} = [];
            Ynew{treeDepth+1} = [];
            [~, TOld_new, ~] = edgeValenceT(TOld, XsToDisconnect);
            [Ynew, Tnew, TRemoved] = YFlipNM_recursive(TOld_new, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect, treeDepth+1);
            % What should I do after I receive the Tnew from other
            % How can I get all the correct new tets?
        end
    end
end
end