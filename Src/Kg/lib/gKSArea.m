function [gs,Ks,Kss]=gKSArea(y1,y2,y3)
% Returns residual and  Jacobian of A=||(y1-y2)x(y2-y3)||
% notation
% YI=Cross_mex(yI)
% q =(y1-y3)x(y2-y3) = Y1y2 - Y1y3 - Y3y2

y1_Crossed = Cross_mex(y1);
y2_Crossed = Cross_mex(y2);
y3_Crossed = Cross_mex(y3);

q= y2_Crossed*y1' - y2_Crossed*y3' + y1_Crossed*y3';
% Q_I = der_yI q =  Cross_mex(yK) -Cross_mex(yJ)
Q1=y2_Crossed-y3_Crossed;
Q2=y3_Crossed-y1_Crossed;
Q3=y1_Crossed-y2_Crossed;
% KK_IJ = der_yJ QI
fact=1/(2*norm(q));
gs=fact.*[Q1'*q; % der_Y1 (det(Y1,Y2,Y3))
    Q2'*q;
    Q3'*q];

Kss=-(2/norm(q)).*(gs)*(gs');

Ks=fact.*[Q1'*Q1               KK(y1_Crossed,y2_Crossed,y3_Crossed,y1,y2,y3) KK(y1_Crossed,y3_Crossed,y2_Crossed,y1,y3,y2);
    KK(y2_Crossed,y1_Crossed,y3_Crossed,y2,y1,y3)  Q2'*Q2              KK(y2_Crossed,y3_Crossed,y1_Crossed,y2,y3,y1);
    KK(y3_Crossed,y1_Crossed,y2_Crossed,y3,y1,y2)  KK(y3_Crossed,y2_Crossed,y1_Crossed,y3,y2,y1) Q3'*Q3            ];
end

