function [Ymat] = Cross(y)
%CROSS Summary of this function goes here
%   Detailed explanation goes here
Ymat = [0 -y(3) y(2)
                y(3) 0 -y(1)
                -y(2) y(1) 0];
end

