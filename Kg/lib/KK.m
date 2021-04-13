function [KIJ]=KK(i,j,k,y1,y2,y3)
Y=[y1;y2;y3];
%KIJ= (Yk-Yj)*(Yi-Yk)+Cross_mex()
KIJ= (Cross_mex(Y(j,:))-Cross_mex(Y(k,:)))*(Cross_mex(Y(i,:))-Cross_mex(Y(k,:)))+...
    Cross_mex((Cross_mex(Y(j,:))*Y(i,:)')')-Cross_mex((Cross_mex(Y(j,:))*Y(k,:)')')-Cross_mex((Cross_mex(Y(k,:))*Y(i,:)')');

end