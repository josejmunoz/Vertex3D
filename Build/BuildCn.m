function [Cn]=BuildCn(T)
%% Building the matrix Cn (Nodal connectivity), bars from tetrahedrons

Cn=zeros(size(T,1)*8,2);
k=0;
for i=1:size(T,1)
     Cn(k+1:k+6,:)=[T(i,1) T(i,2);
                    T(i,1) T(i,3);
                    T(i,1) T(i,4);
                    T(i,2) T(i,3);
                    T(i,2) T(i,4);
                    T(i,3) T(i,4)];
    k=k+6;
end
Cn(k+1:end,:)=[];
end 