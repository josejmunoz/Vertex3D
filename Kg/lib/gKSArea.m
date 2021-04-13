function [gs,Ks,Kss]=gKSArea(y1,y2,y3)
% Returns residual and  Jacobian of A=||(y1-y2)x(y2-y3)||
% notation
% YI=Cross_mex(yI)
% q =(y1-y3)x(y2-y3) = Y1y2 - Y1y3 - Y3y2
q= Cross_mex(y2)*y1' - Cross_mex(y2)*y3' + Cross_mex(y1)*y3';
% Q_I = der_yI q =  Cross_mex(yK) -Cross_mex(yJ)
Q1=Cross_mex(y2)-Cross_mex(y3);
Q2=Cross_mex(y3)-Cross_mex(y1);
Q3=Cross_mex(y1)-Cross_mex(y2);
% KK_IJ = der_yJ QI
fact=1/(2*norm(q));
gs=fact.*[Q1'*q; % der_Y1 (det(Y1,Y2,Y3))
    Q2'*q;
    Q3'*q];

Kss=-(2/norm(q)).*(gs)*(gs');
Ks=fact.*[Q1'*Q1               KK(1,2,3,y1,y2,y3) KK(1,3,2,y1,y2,y3);
    KK(2,1,3,y1,y2,y3)  Q2'*Q2              KK(2,3,1,y1,y2,y3);
    KK(3,1,2,y1,y2,y3)  KK(3,2,1,y1,y2,y3) Q3'*Q3            ];
end

