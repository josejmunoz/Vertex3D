function [contractilityTimeDependent] = contractilityFunction(initC_Contractility, endC_Contractility, initT_Contractility, endT_Contractility, Set)
%CONTRACTILITYFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    slope1 = (endC_Contractility - initC_Contractility) / (endT_Contractility - initT_Contractility);
    contractilityTimeDependent = @(timePoint) (slope1 * (timePoint - Set.TAblation - initT_Contractility)); %initC_Contractility;
end

