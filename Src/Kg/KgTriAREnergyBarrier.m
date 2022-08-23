function [g,K,EnergyB]=KgTriAREnergyBarrier(Geo,Set)
%%KGTRIARENERGYBARRIER Penalise bad aspect ratios
% The residual g and Jacobian K of  Energy Barrier
% Energy  WB =
    [g, K] = initializeKg(Geo, Set);
    EnergyB = 0;
    
    for c=1:Geo.nCells
        if Geo.Remodelling
            if ~ismember(c,Geo.AssembleNodes)
                continue
            end
        end

        if isempty(Geo.Cells(c).AliveStatus) || Geo.Cells(c).AliveStatus ~= 1
            continue
        end

        Cell = Geo.Cells(c);
        Ys = Cell.Y;
        for f = 1:length(Cell.Faces)
            Face = Cell.Faces(f);
            if ~isequal(Face.InterfaceType, 'CellCell')
                Tris = Cell.Faces(f).Tris;
                for t = 1:length(Tris)
                    n3 = Cell.Faces(f).globalIds;
                    nY_original = [Cell.globalIds(Tris(t).Edge)', n3];
                    if Geo.Remodelling
                        if ~any(ismember(nY_original, Geo.AssemblegIds))
                            continue
                        end
                    end

                    y1 = Ys(Tris(t).Edge(1),:);
                    y2 = Ys(Tris(t).Edge(2),:);
                    y3 = Cell.Faces(f).Centre;

                    ys(1, :) = {y1, y2, y3};
                    ys(2, :) = {y2, y3, y1};
                    ys(3, :) = {y3, y1, y2};
                    
                    nY(1, 1:3) = nY_original;
                    nY(2, 1:3) = nY_original([2 3 1]);
                    nY(2, 1:3) = nY_original([3 1 2]);
                    
                    w_t = zeros(3, 1);
                    for numY = 1:size(ys, 1)
                        y1 = ys{numY, 1}';
                        y2 = ys{numY, 2}';
                        y3 = ys{numY, 3}';
                        
                        v_y1 = y2 - y1;
                        v_y2 = y3 - y1;
                        
                        v_y3_1 = y3 - y2;
                        v_y3_2 = y2 - y1;
                        v_y3_3 = -(y3 - y1);

                        w_t(numY) = norm(v_y1)^2 - norm(v_y2)^2;
                        
                        %% g
                        gs(1:3, 1) = Set.lambdaB * w_t(numY)^2 * v_y3_1;
                        gs(4:6, 1) = Set.lambdaB * w_t(numY)^2 * v_y3_2;
                        gs(7:9, 1) = Set.lambdaB * w_t(numY)^2 * v_y3_3;

                        g=Assembleg(g,gs,nY);

                        %% K
                        matrixK = [zeros(3, 3), -eye(3, 3), eye(3, 3);
                            -eye(3, 3), eye(3, 3), zeros(3, 3);
                            eye(3, 3), zeros(3, 3), -eye(3, 3)];
                        
                        Ks = Set.lambdaB * w_t(numY)^2 * matrixK + Set.lambdaB * (gs * gs');

                        K= AssembleK(K,Ks,nY);
                    end
                    EnergyB=EnergyB + Set.lambdaB/2 * sum(w_t);
                end
            end
        end
    end
end
