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
            Tris = Cell.Faces(f).Tris;
            for t = 1:length(Tris)
                n3 = Cell.Faces(f).globalIds;
				nY = [Cell.globalIds(Tris(t).Edge)', n3];
                if Geo.Remodelling
                    if ~any(ismember(nY, Geo.AssemblegIds))
                        continue
                    end
                end
                
                y1 = Ys(Tris(t).Edge(1),:);
                y2 = Ys(Tris(t).Edge(2),:);
                y3 = Cell.Faces(f).Centre;
                
                ys(1, :) = {y1, y2, y3};
                ys(2, :) = {y2, y3, y1};
                ys(3, :) = {y3, y1, y2};
                
                for numY = 1:size(ys, 1)
                    y1 = ys{numY, 1}';
                    y2 = ys{numY, 2}';
                    y3 = ys{numY, 3}';
                    [angle, cos_angle] = ComputeEdgesAngle(y1, y2, y3);
                    cos_angles(numY) = cos_angle;

                    v_y1 = y2 - y1;
                    v_y2 = y3 - y1;

                    w_t = cos_angle - 1/2;
                    beta = 1/(norm(v_y1) * norm(v_y2));

                    gs(1:3, 1) = Set.lambdaB * w_t * (-beta * (v_y1 + v_y2) + cos_angle * (v_y1/(norm(v_y1)^2) + v_y2/(norm(v_y2)^2)));
                    gs(4:6, 1) = Set.lambdaB * w_t * (-beta * (v_y2) + cos_angle * (v_y1/(norm(v_y1)^2)));
                    gs(7:9, 1) = Set.lambdaB * w_t * (-beta * (v_y1) + cos_angle * (v_y2/(norm(v_y2)^2)));

                    g=Assembleg(g,gs,nY);
                    
                    %% K_a
                    K_a(1:3, 1:3) = Set.lambdaB * [2*eye(3) - (v_y1 + v_y2) * (v_y1/(norm(v_y1)^2) + v_y2/(norm(v_y2)^2))'];
                    K_a(4:6, 1:3) = Set.lambdaB * [-eye(3) + v_y2 * (v_y2/(norm(v_y2)^2))'];
                    K_a(7:9, 1:3) = K_a(2, 1);
                    
                    K_a(1:3, 4:6) = Set.lambdaB * [-eye(3) - (v_y1 + v_y2) * (v_y1/(norm(v_y1)^2))'];
                    K_a(4:6, 4:6) = Set.lambdaB * [-v_y2 * (v_y1/(norm(v_y1)^2))'];
                    K_a(7:9, 4:6) = Set.lambdaB * [eye(3) - v_y1 * (v_y1/(norm(v_y1)^2))']; 
                    
                    K_a(1:3, 7:9) = Set.lambdaB * [-eye(3) - (v_y1 + v_y2) * (v_y2/(norm(v_y2)^2))'];
                    K_a(4:6, 7:9) = Set.lambdaB * [-v_y2 * (v_y2/(norm(v_y2)^2))'];
                    K_a(7:9, 7:9) = Set.lambdaB * [eye(3) - v_y1 * (v_y2/(norm(v_y2)^2))'];
                    
                    %% K_b
                    K_b_left(1:3, 1) = v_y1/(norm(v_y1)^2) + v_y2/(norm(v_y2)^2);
                    K_b_left(4:6, 1) = -v_y1/(norm(v_y1)^2);
                    K_b_left(7:9, 1) = -v_y2/(norm(v_y2)^2);
                    
                    K_b_right(1:3, 1) = -Set.lambdaB * (v_y1 + v_y2) + cos_angle * (v_y1/(norm(v_y1)^2) + v_y2/(norm(v_y2)^2));
                    K_b_right(4:6, 1) = -Set.lambdaB * v_y2 + cos_angle * (v_y1/(norm(v_y1)^2));
                    K_b_right(7:9, 1) = -Set.lambdaB * v_y2 + cos_angle * (v_y2/(norm(v_y2)^2));
                    
                    K_b = K_b_left * (K_b_right)';
                    
                    %% K_c
                    K_c(1:3, 1:3) = -(1/norm(v_y1)^2 + 1/(norm(v_y2)^2))*eye(3) + (v_y1 * v_y1')/norm(v_y1)^4 - (v_y2 * v_y2')/norm(v_y2)^4;
                    K_c(4:6, 1:3) = (1/norm(v_y1)^2)*eye(3) - (v_y1 * v_y1')/norm(v_y1)^4;
                    K_c(7:9, 1:3) = (1/norm(v_y2)^2)*eye(3) - (v_y2 * v_y2')/norm(v_y2)^4;
                    
                    K_c(1:3, 4:6) = (1/norm(v_y1)^2)*eye(3) - (v_y1 * v_y1')/norm(v_y1)^4;
                    K_c(4:6, 4:6) = -(1/norm(v_y1)^2)*eye(3) + (v_y1 * v_y1')/norm(v_y1)^4;
                    K_c(7:9, 4:6) = 0;
                    
                    K_c(1:3, 7:9) = (1/norm(v_y2)^2)*eye(3) - (v_y2 * v_y2')/norm(v_y2)^4;
                    K_c(4:6, 7:9) = 0;
                    K_c(7:9, 7:9) = -(1/norm(v_y2)^2)*eye(3) + (v_y2 * v_y2')/norm(v_y2)^4;
                    
                    %% Ks
                    Ks = Set.lambdaB * (gs * gs') + Set.lambdaB * w_t * (K_a + K_b + cos_angle * K_c);

                    K= AssembleK(K,Ks,nY);
                end
                EnergyB=EnergyB + Set.lambdaB/2 * sum((cos_angles - 1/2).^2);
            end
        end
    end
end
