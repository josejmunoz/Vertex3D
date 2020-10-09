function [g,Cell,EnergyBend]=gBending(Cell,Y,Set)
% K(i,j)= derivative of g(i) wrt to x(j)
% energy based on the total cell area W_s= sum_cell ((As-As0)/As0)^2

%% Input 
Set.Sparse=true;

%% Set parameters
ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual



EnergyBend=0;

%% Compute Volume
% [Cell]=ComputeCellVolume(Cell,Y);
for i=1:ncell
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes) 
           continue
        end 
    end 
    L=Set.lambdaBend;
    Edges=Cell.Edges{i};
    for e=1:size(Edges,1)
        if Edges(e,1)<=Y.n, Y1=Y.DataRow(Edges(e,1),:); else,  Y1=Cell.SurfsCenters.DataRow(Edges(e,1)-Y.n,:); end 
        if Edges(e,2)<=Y.n, Y2=Y.DataRow(Edges(e,2),:); else,  Y2=Cell.SurfsCenters.DataRow(Edges(e,2)-Y.n,:); end 
        if Edges(e,3)<=Y.n, Y3=Y.DataRow(Edges(e,3),:); else,  Y3=Cell.SurfsCenters.DataRow(Edges(e,3)-Y.n,:); end 
        if Edges(e,4)<=Y.n, Y4=Y.DataRow(Edges(e,4),:); else,  Y4=Cell.SurfsCenters.DataRow(Edges(e,4)-Y.n,:);end 
        
        n1=cross(Y2-Y1,Y3-Y1);  n2=cross(Y4-Y1,Y2-Y1);
        A1=(1/2)*norm(n1);      A2=(1/2)*norm(n2);
        B=dot(n1,n2);
        n1=n1/norm(n1);        n2=n2/norm(n2);
        nn=dot(n1,n2);
        if nn==-1
            nn=1;
        end 
        
        fact0=(L/12)*(1-(2/sqrt(2*nn+2)));
        fact1=(L/6)*(1-sqrt(2*nn+2)/2)^2;
        
        [dA1dy,dA2dy]=dAdY(Y1,Y2,Y3,Y4);
        [dBdy]=dBdY(Y1,Y2,Y3,Y4);
        dAdy=dA1dy+dA2dy;
        
        S=( (1/B)*dBdy - (1/A1)*dA1dy -(1/A2)*dA2dy );
        dnndy=nn*S;
        
        if Set.BendingAreaDependent       
            ge=fact0*(A1+A2)*dnndy + fact1*dAdy;
        else
            ge=fact0*dnndy;
        end 
        % Assembleg
        dim=3;
        for I=1:length(Edges(e,:)) % loop on 3 vertices of triangle
                idofg=(Edges(e,I)-1)*dim+1:Edges(e,I)*dim; % global dof
                idofl=(I-1)*dim+1:I*dim;
                g(idofg)=g(idofg)+ge(idofl);
        end
        

    end 




end 





end
%%
function [dBdy]=dBdY(y1,y2,y3,y4)
y12=y2'-y1'; y23=y3'-y2'; y24=y4'-y2'; y13=y3'-y1'; y14=y4'-y1'; 
dBdy1= Cross(y24)*Cross(y12)*y13-Cross(y23)*Cross(y14)*y12;
dBdy2=-Cross(y14)*Cross(y12)*y13+Cross(y13)*Cross(y14)*y12;
dBdy3=-Cross(y12)*Cross(y14)*y12;
dBdy4= Cross(y12)*Cross(y12)*y13;
dBdy=[dBdy1;dBdy2;dBdy3;dBdy4]; 





end 

function [ga1,ga2]=dAdY(y1,y2,y3,y4)

q1= Cross(y2)*y1' - Cross(y2)*y3' + Cross(y1)*y3';
q2= Cross(y2)*y1' - Cross(y2)*y4' + Cross(y1)*y4';

Q11=Cross(y2)-Cross(y3);
Q21=Cross(y3)-Cross(y1);
Q12=Cross(y2)-Cross(y4);
Q22=Cross(y4)-Cross(y1);
Q3=Cross(y1)-Cross(y2);
fact1=1/(2*norm(q1));
fact2=1/(2*norm(q2));

ga1=[fact1*Q11'*q1;
     fact1*Q21'*q1;
     fact1*Q3'*q1;
     zeros(3,1)];
ga2=[fact2*Q12'*q2;
     fact2*Q22'*q2;
     zeros(3,1);
     fact2*Q3'*q2];

end

%%

function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];

end

