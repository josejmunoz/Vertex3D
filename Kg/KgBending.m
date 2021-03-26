function [g,K,Cell,Energy]=KgBending(Cell,Y,Set)
% The residual g and Jacobian K of the Bending Energy
%  Potential: -> Set.BendingAreaDependent=1 : Wb=(1/2) lambdaBend* sum_edges( 1-cos(theta/2)^2*(At1+At2)
%             -> Set.BendingAreaDependent=0 : Wb=(1/2) lambdaBend* sum_edges( 1-cos(theta/2)^2
%                       where  theta: the angle between the pair of triangles
%                              At1 and At2 : the area of the triangles

ncell=Cell.n;

%% Initialize
dimg=Set.NumTotalV*3;

g=zeros(dimg,1); % Local cell residual
if Set.Sparse && nargout>1
    sk=0;
    si=zeros((dimg*3)^2,1); % Each vertex is shared by at least 3 cells
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
elseif nargout>1
    K=zeros(dimg); % Also used in sparse
end


Energy=0;


%% Loop over cells 
for i=1:ncell
    if Cell.DebrisCells(i)
        continue;
    end 
    if ~Cell.AssembleAll
        if ~ismember(Cell.Int(i),Cell.AssembleNodes)
            continue
        end
    end
    L=Set.lambdaBend;
    Edges=Cell.Edges{i};
    for e=1:size(Edges,1)
        if ~Cell.AssembleAll && ~any(ismember(Edges(e,:),Cell.RemodelledVertices))
            continue
        end
        if Edges(e,1)<=Y.n, Y1=Y.DataRow(Edges(e,1),:); else,  Y1=Cell.FaceCentres.DataRow(Edges(e,1)-Y.n,:); end
        if Edges(e,2)<=Y.n, Y2=Y.DataRow(Edges(e,2),:); else,  Y2=Cell.FaceCentres.DataRow(Edges(e,2)-Y.n,:); end
        if Edges(e,3)<=Y.n, Y3=Y.DataRow(Edges(e,3),:); else,  Y3=Cell.FaceCentres.DataRow(Edges(e,3)-Y.n,:); end
        if Edges(e,4)<=Y.n, Y4=Y.DataRow(Edges(e,4),:); else,  Y4=Cell.FaceCentres.DataRow(Edges(e,4)-Y.n,:);end
        
        
        
        
        
        n1=cross(Y2-Y1,Y3-Y1);  n2=cross(Y4-Y1,Y2-Y1);
        A1=(1/2)*norm(n1);      A2=(1/2)*norm(n2);
        B=dot(n1,n2);
        n1=n1/norm(n1);        n2=n2/norm(n2);
        nn=dot(n1,n2);
        if nn==-1
            nn=1;
        end
        fact0=(L/12)*(1-2/sqrt(2*nn+2));
        fact1=(L/6)*(1-sqrt(2*nn+2)/2)^2;
        fact2=(L/6)*(1 / ((2*nn+2)^(3/2)) );
        
        [dA1dy,dA2dy,dA1dydy,dA2dydy]=dAdY(Y1,Y2,Y3,Y4);
        [dBdy,dBdydy]=dBdY(Y1,Y2,Y3,Y4);
        dAdy=dA1dy+dA2dy;
        
        
        if abs(B)<eps
            S=zeros(size(dBdy));
            dnndy=zeros(size(S));
            Knn=zeros(length(S));
        else
            S=( (1/B)*dBdy - (1/A1)*dA1dy -(1/A2)*dA2dy );
            dnndy=nn*S;
            Knn=nn * ( S*S' + (1/B)*dBdydy - (1/B^2)*dBdy*(dBdy')...
                -(1/A1)*dA1dydy +(1/A1^2)*dA1dy*(dA1dy')...
                -(1/A2)*dA2dydy +(1/A2^2)*dA2dy*(dA2dy'));
        end
        
        if Set.BendingAreaDependent
            ge=fact0*(A1+A2)*dnndy + fact1*dAdy;
            Ke=fact2*(A1+A2)*dnndy*(dnndy') +  fact0*( (A1+A2)*Knn + dnndy*dAdy' + dAdy*dnndy') + fact1*(dA1dydy+dA2dydy);
        else
            % Area independent
            ge=fact0*dnndy;
            Ke=fact2*dnndy*(dnndy')+  fact0*Knn;
        end
        % Assembleg
        dim=3;
        for I=1:length(Edges(e,:)) % loop on 3 vertices of triangle
            idofg=(Edges(e,I)-1)*dim+1:Edges(e,I)*dim; % global dof
            idofl=(I-1)*dim+1:I*dim;
            g(idofg)=g(idofg)+ge(idofl);
        end
        
        % AssembleK
        if  nargout>1
            if Set.Sparse
                [si,sj,sv,sk]= AssembleKSparse(Ke,Edges(e,:),si,sj,sv,sk);
            else
                K= AssembleK(K,Ke,Edges(e,:));
            end
            Energy=Energy+(L/6)*(1-sqrt(2*nn+2)/2)^2;
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg);
end






%% Loop over Cells
%     % Analytical residual g and Jacobian K
% for i=1:ncell
%     if ~Cell.AssembleAll
%         if ~ismember(Cell.Int(i),Cell.AssembleNodes)
%            continue
%         end
%     end
%     lambdaS=Set.lambdaS;
% %     fact=( lambdaS / (4*Cell.SArea(i)) ) *(  (Cell.SArea(i)-Cell.SArea0(i)) / Cell.SArea0(i)^2   );
% if Set.A0eq0
%     fact=lambdaS *  (Cell.SArea(i)) / Cell.SArea0(i)^2   ;
% else
%     fact=lambdaS *  (Cell.SArea(i)-Cell.SArea0(i)) / Cell.SArea0(i)^2   ;
% end
%     ge=zeros(dimg,1); % Local cell residual
% %             K2=zeros(dimg); % Also used in sparse
%
%     % Loop over Cell-face-triangles
%     Tris=Cell.Tris{i};
%     for t=1:size(Tris,1)
%         nY=Tris(t,:);
%         Y1=Y.DataRow(nY(1),:);
%         Y2=Y.DataRow(nY(2),:);
%         Y4=Cell.FaceCentres.DataRow(nY(3),:);
%         nY(3)=nY(3)+Set.NumMainV;
%         [gs,Ks,Kss]=gKSArea(Y1,Y2,Y4);
%         Ks=fact*(Ks+Kss);
%         ge=AssemblegSArea(ge,gs,nY);
%         if Set.Sparse
%             [si,sj,sv,sk]= AssembleKSparse(Ks,nY,si,sj,sv,sk);
%
%         else
%             K= AssembleK(K,Ks,nY);
%         end
%     end
%
%     g=g+ge*fact; % Volume contribution of each triangle is (y1-y2)'*J*(y2-y3)/2
%     if Set.Sparse
%         K=K+sparse((ge)*(ge')*lambdaS/(Cell.SArea0(i)^2));
%     else
%        K=K+(ge)*(ge')*lambdaS/(Cell.SArea0(i)^2);  %-(gee)*(gee')*(fact);
%     end
%
%     if Set.A0eq0
%         EnergyBend=EnergyBend+ lambdaS/2 *((Cell.SArea(i)) / Cell.SArea0(i))^2;
%     else
%         EnergyBend=EnergyBend+ lambdaS/2 *((Cell.SArea(i)-Cell.SArea0(i)) / Cell.SArea0(i))^2;
%     end
%
% end
%
% if Set.Sparse
%     K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
% end
end
%%
function [dBdy,dBdydy]=dBdY(y1,y2,y3,y4)
y12=y2'-y1'; y23=y3'-y2'; y24=y4'-y2'; y13=y3'-y1'; y14=y4'-y1';
dBdy1= Cross(y24)*Cross(y12)*y13-Cross(y23)*Cross(y14)*y12;
dBdy2=-Cross(y14)*Cross(y12)*y13+Cross(y13)*Cross(y14)*y12;
dBdy3=-Cross(y12)*Cross(y14)*y12;
dBdy4= Cross(y12)*Cross(y12)*y13;
dBdy=[dBdy1;dBdy2;dBdy3;dBdy4];

dBdy1dy1=Cross(y24)*Cross(y23)+Cross(y23)*Cross(y24);
dBdy2dy2=Cross(y14)*Cross(y13)+Cross(y13)*Cross(y14);
dBdy3dy3=zeros(3);
dBdy4dy4=zeros(3);

dBdy1dy2=-Cross(y24)*Cross(y13) + Cross(Cross(y12)*y13) - Cross(y23)*Cross(y14) - Cross(Cross(y14)*y12);
dBdy1dy3= Cross(y24)*Cross(y12) + Cross(Cross(y14)*y12);
dBdy1dy4= Cross(y23)*Cross(y12) - Cross(Cross(y12)*y13);

dBdy2dy3=-Cross(y14)*Cross(y12) - Cross(Cross(y14)*y12);
dBdy2dy4=-Cross(y13)*Cross(y12) + Cross(Cross(y12)*y13);

dBdy3dy4=Cross(y12)*Cross(y12);


dBdydy=[dBdy1dy1  dBdy1dy2  dBdy1dy3  dBdy1dy4;
    dBdy1dy2' dBdy2dy2  dBdy2dy3  dBdy2dy4;
    dBdy1dy3' dBdy2dy3' dBdy3dy3  dBdy3dy4;
    dBdy1dy4' dBdy2dy4' dBdy3dy4' dBdy4dy4];
end

function [ga1,ga2,Ka1,Ka2]=dAdY(y1,y2,y3,y4)

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
Ka1=-(2/norm(q1)).*(ga1)*(ga1');
Ka2=-(2/norm(q2)).*(ga2)*(ga2');
Kaa1=fact1.*[Q11'*Q11               KK(1,2,3,y1,y2,y3)    KK(1,3,2,y1,y2,y3) zeros(3);
    KK(2,1,3,y1,y2,y3)   Q21'*Q21                KK(2,3,1,y1,y2,y3) zeros(3);
    KK(3,1,2,y1,y2,y3)   KK(3,2,1,y1,y2,y3)    Q3'*Q3             zeros(3);
    zeros(3)             zeros(3)              zeros(3)           zeros(3)];

Kaa2=fact2.*[Q12'*Q12               KK(1,2,3,y1,y2,y4)    zeros(3)    KK(1,3,2,y1,y2,y4);
    KK(2,1,3,y1,y2,y4)   Q22'*Q22                zeros(3)    KK(2,3,1,y1,y2,y4);
    zeros(3)             zeros(3)              zeros(3)    zeros(3)
    KK(3,1,2,y1,y2,y4)   KK(3,2,1,y1,y2,y4)    zeros(3)    Q3'*Q3];
Ka1=Ka1+Kaa1;
Ka2=Ka2+Kaa2;
end

