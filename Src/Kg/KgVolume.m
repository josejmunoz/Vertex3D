function [g,K,EnergyV]=KgVolume(Geo, Set)
	% The residual g and Jacobian K of Volume Energy 
	% Energy W_s= sum_cell lambdaV ((V-V0)/V0)^2
	[g, K] = initializeKg(Geo, Set);
	EnergyV = 0;
	
	%% Loop over Cells 
	% Analytical residual g and Jacobian K
	for c=1:Geo.nCells
		if Geo.Remodelling
			if ~ismember(c,Geo.AssembleNodes)
        		continue
			end
		end
		Cell = Geo.Cells(c);
		Ys = Cell.Y;
    	lambdaV=Set.lambdaV;
    	fact=lambdaV*(Cell.Vol-Cell.Vol0)/Cell.Vol0^2;
    	
    	ge=sparse(size(g, 1), 1);
		for f = 1:length(Cell.Faces)
			Tris = Cell.Faces(f).Tris;
			for t=1:length(Tris)
				y1 = Ys(Tris(t,1),:);
				y2 = Ys(Tris(t,2),:);
				if length(Tris) == 3
					y3 = Ys(Tris(t+1,2),:);
					n3 = Cell.globalIds(Tris(t+1,2));
				else
					y3 = Cell.Faces(f).Centre;
					n3 = Cell.Faces(f).globalIds;
				end
				nY = [Cell.globalIds(Tris(t,:))', n3];
				if Geo.Remodelling
					if ~any(ismember(nY,Geo.AssemblegIds))
                        if length(Tris) == 3
                            break
                        else
                		    continue
                        end
					end
				end
				[gs,Ks]=gKDet(y1, y2, y3); 
				ge=Assembleg(ge,gs,nY); 
% 				fprintf("%.15f %.2f %.2f %.2f %d %d %d \n", norm(g), y3, c, f, t);
				K = AssembleK(K,Ks*fact/6,nY);
				if length(Tris) == 3
					break
				end
			end
		end
    	g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
    	geMatrix = lambdaV*((ge)*(ge')/6/6/Cell.Vol0^2);
    	K=K+geMatrix;
    	EnergyV=EnergyV+lambdaV/2 *((Cell.Vol-Cell.Vol0)/Cell.Vol0)^2;    
    end
end
%%
% TODO FIXME: Move this to lib?
function [gs,Ks]=gKDet(Y1,Y2,Y3)
	% Returns residual and  Jacobian of det(Y)=y1'*cross(y2,y3)
	% gs=[der_y1 det(Y) der_y2 det(Y) der_y3 det(Y)]
	% Ks=[der_y1y1 det(Y) der_y1y2 det(Y) der_y1y3 det(Y)
	%     der_y2y1 det(Y) der_y2y2 det(Y) der_y2y3 det(Y)
	%     der_y3y1 det(Y) der_y3y2 det(Y) der_y3y3 det(Y)]
	dim=length(Y1);
	gs=[cross(Y2,Y3)'; % der_Y1 (det(Y1,Y2,Y3)) 
    	cross(Y3,Y1)';
    	cross(Y1,Y2)'];
	Ks=[ zeros(dim) -Cross_mex(Y3)   Cross_mex(Y2) % g associated to der wrt vertex 1
    	Cross_mex(Y3)   zeros(dim) -Cross_mex(Y1)
    	-Cross_mex(Y2)   Cross_mex(Y1)  zeros(dim)];
end


