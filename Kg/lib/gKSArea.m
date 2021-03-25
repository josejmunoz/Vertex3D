function [gs,Ks,Kss] = gKSArea(y1,y2,y3)
%GSKAREA Returns residual and  Jacobian of A=||(y1-y2)x(y2-y3)||
%   notation
%   YI=Cross(yI)
%   q =(y1-y3)x(y2-y3) = Y1y2 - Y1y3 - Y3y2

y2_Cross = Cross(y2);
y1_Cross = Cross(y1);
y3_Cross = Cross(y3);
q = y2_Cross*y1' - y2_Cross*y3' + y1_Cross*y3';

% Q_I = der_yI q =  Cross(yK) -Cross(yJ)
Q1 = y2_Cross - y3_Cross;
Q2 = y3_Cross - y1_Cross;
Q3 = y1_Cross - y2_Cross;

% KK_IJ = der_yJ QI
fact=1/(2*norm(q));
gs=fact.*[Q1'*q; % der_Y1 (det(Y1,Y2,Y3))
    Q2'*q;
    Q3'*q];

Kss=-(2/norm(q)).*(gs)*(gs');

Y = [y1;y2;y3];

Ks=fact.*[Q1'*Q1 KK(1,2,3,Y) KK(1,3,2,Y);
    KK(2,1,3,Y) Q2'*Q2 KK(2,3,1,Y);
    KK(3,1,2,Y) KK(3,2,1,Y) Q3'*Q3];
end

