function [obj]=ComputeTriArea(obj,Y,SurfsCenters)
    for c = 1:3
        Cell = Geo.Cells(c);
        Ys = Cell.Y;
        for f = 1:length(Cell.Faces)
            Face = Cell.Faces(f);
            for t = 1:length(Face.Tris)
                Tri = Face.Tris(t);
                YTri = [Ys(Tri(1),:), Ys(Tri(2),:), Face.Centre];
                area = (1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                nrg  = exp(Set.lambdaB*(1-Set.Beta*area/Set.BarrierTri0));
            end
        end
    end

    for i=1:obj.n
        if obj.NotEmpty(i)
            if length(obj.Vertices{i})==3
                YTri=Y(obj.Vertices{i},:);
                obj.AreaTri{i}=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                %                      obj.AreaTri0{i}=1e-3;
                obj.Area(i)=sum(obj.AreaTri{i});
            else 
                obj.AreaTri{i}=zeros(length(obj.Vertices{i}),1);
            %                        obj.AreaTri0{i}=ones(length(obj.Vertices{i}),1)*1e-3;
                for j=1:length(obj.Vertices{i})-1
                    YTri=[Y(obj.Vertices{i}([j j+1]),:); SurfsCenters(i,:)];
                    obj.AreaTri{i}(j)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                end
                YTri=[Y(obj.Vertices{i}([j+1 1]),:); SurfsCenters(i,:)];
                obj.AreaTri{i}(j+1)=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
                obj.Area(i)=sum(obj.AreaTri{i});
            end 
        else 
            obj.AreaTri{i}=[];
        %                    obj.AreaTri0{i}=[];
        end
    end
end 