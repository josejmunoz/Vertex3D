function [Geo, Set] = readBatchLine(tlines, numLine, Set, Geo)
%READBATCHLINE Summary of this function goes here
%   Detailed explanation goes here
eval(tlines{numLine});
end

