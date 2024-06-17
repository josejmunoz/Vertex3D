function [circularity] = Compute2DCircularity(area, perimeter)
%COMPUTE2DCIRCULARITY Calculate circularity
% 
circularity = (4 * pi * area) / (perimeter^2);
end

