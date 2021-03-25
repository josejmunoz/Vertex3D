function [KIJ]=KK(i,j,k, Y)
%KK Summary of this function goes here
%   Detailed explanation goes here

%KIJ= (Yk-Yj)*(Yi-Yk)+Cross()
Y_cross_i = Cross(Y(i,:));
Y_cross_j = Cross(Y(j,:));
Y_cross_k = Cross(Y(k,:));

KIJ= (Y_cross_j - Y_cross_k) * (Y_cross_i - Y_cross_k) + ...
    Cross(Y_cross_j*Y(i,:)')-Cross(Y_cross_j * Y(k,:)')- Cross(Y_cross_k * Y(i,:)');

end

