function finalValue = AddNoiseToParameter(avgParameter, noise, currentTri)
minValue = avgParameter - avgParameter*noise;
maxValue = avgParameter + avgParameter*noise;
if minValue < 0
    minValue = eps;
end
finalValue = minValue + (maxValue-minValue)*rand();

% if exist('currentTri', 'var')
%     if ~isempty(currentTri.pastContractilityValue)
%         finalValue = mean([finalValue, currentTri.pastContractilityValue]);
%     end
% end
end