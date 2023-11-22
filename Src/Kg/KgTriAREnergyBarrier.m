function [g,K,Energy_T]=KgTriAREnergyBarrier(Geo,Set)
%%KGTRIARENERGYBARRIER Penalise bad aspect ratios
% The residual g and Jacobian K of  Energy Barrier
% Energy  WB =
    [g, K] = initializeKg(Geo, Set);
    IDs=[Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    Energy=zeros(max(IDs),1);
    fact=Set.lambdaR/Set.lmin0^4;
    for c = IDs
        if Geo.Remodelling
            if ~ismember(c,Geo.AssembleNodes)
                continue
            end
        end

        if Geo.Cells(c).AliveStatus
            Cell = Geo.Cells(c);
            Ys = Cell.Y;
            for f = 1:length(Cell.Faces)
                Face = Cell.Faces(f);
                if ~isequal(Face.InterfaceType, 2)
                    Tris = Cell.Faces(f).Tris;
                    for t = 1:length(Tris)
                        n3 = Cell.Faces(f).globalIds;
                        nY_original = [Cell.globalIds(Tris(t).Edge)', n3];
                        if Geo.Remodelling
                            if ~any(ismember(nY_original, Geo.AssemblegIds))
                                continue
                            end
                        end
    
                        y1 = Ys(Tris(t).Edge(1),:)';
                        y2 = Ys(Tris(t).Edge(2),:)';
                        y3 = Cell.Faces(f).Centre';
    
                        y12=y1-y2;
                        y23=y2-y3;
                        y31=y3-y1;
                        
                        w1=norm(y31)^2-norm(y12)^2;
                        w2=norm(y12)^2-norm(y23)^2;
                        w3=norm(y23)^2-norm(y31)^2;
                        
                        g1=[y23;y12;y31];
                        g2=[y12;y31;y23];
                        g3=[y31;y23;y12];
                        
                        gs=2*(w1*g1+w2*g2+w3*g3);
                        
                        I=eye(3);
                        Ks=2*[(w2-w3)*I (w1-w2)*I (w3-w1)*I
                              (w1-w2)*I (w3-w1)*I (w2-w3)*I
                              (w3-w1)*I (w2-w3)*I (w1-w2)*I]+4*(g1*g1'+g2*g2'+g3*g3');
                                                    
                        g=Assembleg(g,gs*fact,nY_original);
                        K= AssembleK(K,Ks*fact,nY_original);
                        
                        Energy(c)=Energy(c)+fact/2*(w1^2+w2^2+w3^2);
                    end
                end
            end
        end
    end
    Energy_T = sum(Energy);
end
