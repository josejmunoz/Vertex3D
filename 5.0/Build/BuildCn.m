function [Cn]=BuildCn(Twg,XgID)


Cn=zeros(size(Twg,1)*8,2);
    k=0;
    
%% bars  from tetrahedrons
for i=1:size(Twg,1)
    if abs(sum(ismember(Twg(i,:),XgID))-4)>eps
     Cn(k+1:k+6,:)=[Twg(i,1) Twg(i,2);
                    Twg(i,1) Twg(i,3);
                    Twg(i,1) Twg(i,4);
                    Twg(i,2) Twg(i,3);
                    Twg(i,2) Twg(i,4);
                    Twg(i,3) Twg(i,4)];
    k=k+6;
    
    end 
end

 Cn(k+1:end,:)=[];

end 