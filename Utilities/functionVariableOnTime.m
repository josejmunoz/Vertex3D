function [contractilityTimeDependent] = functionVariableOnTime(initValue, endValue, initTime, endTime, Set)
%FUNCTIONVARIABLEONTIME Summary of this function goes here
%   Detailed explanation goes here
    slope1 = (endValue - initValue) / (endTime - initTime);
    contractilityTimeDependent = @(timePoint) (slope1 * (timePoint - Set.TAblation - initTime) + initValue);
end

