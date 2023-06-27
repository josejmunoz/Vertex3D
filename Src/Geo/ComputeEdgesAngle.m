function [angle, cos_angle] = ComputeEdgesAngle(y1, y2, y3)
%COMPUTEEDGESANGLE Summary of this function goes here
%   Detailed explanation goes here

v_y1 = y2 - y1;
v_y2 = y3 - y1;

cos_angle = (v_y1' * v_y2) / (norm(v_y1) * norm(v_y2));

angle = acos (cos_angle);

end

