function [Cell,Y]=CheckOrderingOfTriangulaiton(Cell,Y,Set)
%% This function makes sure that the triangulation of cell-Faces is ordered correctly.   


    Sides=[1 2;
           2 3;
           3 1];
    Recompute=false;
    for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        Cell = Geo.Cells(c);
        IsConsistent=true;
        for f = 1:length(Cell.Faces)
            Face = Cell.Faces(f);
            Tris = Face.Tris;
            Tris(Tris(:,3)>0,3)=Tris(Tris(:,3)>0,3)+Geo.numY;
            Tris(Tris(:,3)<0,3)=abs(Tris(Tris(:,3)<0,3));
            Id=1:size(Tris,1);
            for t = 1:length(Tris)
                Side1=Tris(t,Sides(1,:));
                TID=Id(sum(ismember(Tris,Side1),2)==2); TID(TID==t)=[];
                if (Side1(1)==Tris(TID,1) && Side1(2)==Tris(TID,2)) || ...
                   (Side1(1)==Tris(TID,2) && Side1(2)==Tris(TID,3)) || ...
                   (Side1(1)==Tris(TID,3) && Side1(2)==Tris(TID,1))
        %            fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
                   IsConsistent=false;
                   break
                end 
                Side2=Tris(t,Sides(2,:));
                TID=Id(sum(ismember(Tris,Side2),2)==2); TID(TID==t)=[];
                if (Side2(1)==Tris(TID,1) && Side2(2)==Tris(TID,2)) || ...
                   (Side2(1)==Tris(TID,2) && Side2(2)==Tris(TID,3)) || ...
                   (Side2(1)==Tris(TID,3) && Side2(2)==Tris(TID,1))
        %            fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
                   IsConsistent=false;
                   break
                end
                Side3=Tris(t,Sides(3,:));
                TID=Id(sum(ismember(Tris,Side3),2)==2); TID(TID==t)=[];
                if (Side3(1)==Tris(TID,1) && Side3(2)==Tris(TID,2)) || ...
                   (Side3(1)==Tris(TID,2) && Side3(2)==Tris(TID,3)) || ...
                   (Side3(1)==Tris(TID,3) && Side3(2)==Tris(TID,1))
        %            fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
                   IsConsistent=false;
                   break
                end

                if IsConsistent,  continue; end 
                warning('Inconsistent triangulation was detected.... \n')
                Recompute=true;
                    
                % Choose a particular consistent order between the two options 
                Included=zeros(size(Tris,1),1);
                TrisAux=zeros(size(Tris));
                TrisAux(1,:)=Tris(1,:);
                Included(1)=1;
                NotAllIncluded=true;
                while NotAllIncluded
                    for t=2:size(Tris,1)
                        if Included(t)==1, continue; end
                        T=Tris(t,:);
                        for tt=1:size(Tris,1)
                            if tt==t || Included(tt)==0, continue; end 
                            TT=[TrisAux(tt,[1 2]) Tris(tt,3)];
                            if (T(1)==TT(1) && T(2)==TT(2)) || ...
                               (T(1)==TT(2) && T(2)==TT(3)) || ...
                               (T(1)==TT(3) && T(2)==TT(1)) || ...
                               (T(2)==TT(1) && T(3)==TT(2)) || ...
                               (T(2)==TT(2) && T(3)==TT(3)) || ...
                               (T(2)==TT(3) && T(3)==TT(1)) || ...
                               (T(3)==TT(1) && T(1)==TT(2)) || ...
                               (T(3)==TT(2) && T(1)==TT(3)) || ...
                               (T(3)==TT(3) && T(1)==TT(1))
                                TrisAux(t,:)=[Cell.Tris{c}(t,2) Cell.Tris{c}(t,1) Cell.Tris{c}(t,3)];
                                Included(t)=1;
                                break
                            elseif (T(2)==TT(1) && T(1)==TT(2)) || ...
                                   (T(2)==TT(2) && T(1)==TT(3)) || ...
                                   (T(2)==TT(3) && T(1)==TT(1)) || ...
                                   (T(3)==TT(1) && T(2)==TT(2)) || ...
                                   (T(3)==TT(2) && T(2)==TT(3)) || ...
                                   (T(3)==TT(3) && T(2)==TT(1)) || ...
                                   (T(1)==TT(1) && T(3)==TT(2)) || ...
                                   (T(1)==TT(2) && T(3)==TT(3)) || ...
                                   (T(1)==TT(3) && T(3)==TT(1)) 
                                   TrisAux(t,:)=Cell.Tris{c}(t,:);
                                   Included(t)=1;
                                   break
                            end 
                        end 
                    end
                    if sum(Included)==length(Included)
                        TrisAux(1,:)=Cell.Tris{c}(1,:);
                        NotAllIncluded=false;
                    end
                end
            end
        end
        
    end


    % loop over Cell (Checking the order of vertices of =each triangle )
    for c=1:Cell.n 
        IsConsistent=true;
        Tris=Cell.Tris{c};
        Tris(Tris(:,3)>0,3)=Tris(Tris(:,3)>0,3)+Y.n;
        Tris(Tris(:,3)<0,3)=abs(Tris(Tris(:,3)<0,3));
        Id=1:size(Tris,1);
        for t=1:size(Tris,1)
            Side1=Tris(t,Sides(1,:));
            TID=Id(sum(ismember(Tris,Side1),2)==2); TID(TID==t)=[];
            if (Side1(1)==Tris(TID,1) && Side1(2)==Tris(TID,2)) || ...
               (Side1(1)==Tris(TID,2) && Side1(2)==Tris(TID,3)) || ...
               (Side1(1)==Tris(TID,3) && Side1(2)==Tris(TID,1))
    %            fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
               IsConsistent=false;
               break
            end 
            Side2=Tris(t,Sides(2,:));
            TID=Id(sum(ismember(Tris,Side2),2)==2); TID(TID==t)=[];
            if (Side2(1)==Tris(TID,1) && Side2(2)==Tris(TID,2)) || ...
               (Side2(1)==Tris(TID,2) && Side2(2)==Tris(TID,3)) || ...
               (Side2(1)==Tris(TID,3) && Side2(2)==Tris(TID,1))
    %            fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
               IsConsistent=false;
               break
            end
            Side3=Tris(t,Sides(3,:));
            TID=Id(sum(ismember(Tris,Side3),2)==2); TID(TID==t)=[];
            if (Side3(1)==Tris(TID,1) && Side3(2)==Tris(TID,2)) || ...
               (Side3(1)==Tris(TID,2) && Side3(2)==Tris(TID,3)) || ...
               (Side3(1)==Tris(TID,3) && Side3(2)==Tris(TID,1))
    %            fprintf('Triangles %i and %i of %i cell are incompatible \n',t,TID,c);
               IsConsistent=false;
               break
            end
        end 
        
        if IsConsistent,  continue; end 
        warning('Inconsistent triangulation was detected.... \n')
        Recompute=true;
            
        % Choose a particular consistent order between the two options 
        Included=zeros(size(Tris,1),1);
        TrisAux=zeros(size(Tris));
        TrisAux(1,:)=Tris(1,:);
        Included(1)=1;
        NotAllIncluded=true;
        while NotAllIncluded
            for t=2:size(Tris,1)
                if Included(t)==1, continue; end
                T=Tris(t,:);
                for tt=1:size(Tris,1)
                    if tt==t || Included(tt)==0, continue; end 
                    TT=[TrisAux(tt,[1 2]) Tris(tt,3)];
                    if (T(1)==TT(1) && T(2)==TT(2)) || ...
                       (T(1)==TT(2) && T(2)==TT(3)) || ...
                       (T(1)==TT(3) && T(2)==TT(1)) || ...
                       (T(2)==TT(1) && T(3)==TT(2)) || ...
                       (T(2)==TT(2) && T(3)==TT(3)) || ...
                       (T(2)==TT(3) && T(3)==TT(1)) || ...
                       (T(3)==TT(1) && T(1)==TT(2)) || ...
                       (T(3)==TT(2) && T(1)==TT(3)) || ...
                       (T(3)==TT(3) && T(1)==TT(1))
                        TrisAux(t,:)=[Cell.Tris{c}(t,2) Cell.Tris{c}(t,1) Cell.Tris{c}(t,3)];
                        Included(t)=1;
                        break
                    elseif (T(2)==TT(1) && T(1)==TT(2)) || ...
                           (T(2)==TT(2) && T(1)==TT(3)) || ...
                           (T(2)==TT(3) && T(1)==TT(1)) || ...
                           (T(3)==TT(1) && T(2)==TT(2)) || ...
                           (T(3)==TT(2) && T(2)==TT(3)) || ...
                           (T(3)==TT(3) && T(2)==TT(1)) || ...
                           (T(1)==TT(1) && T(3)==TT(2)) || ...
                           (T(1)==TT(2) && T(3)==TT(3)) || ...
                           (T(1)==TT(3) && T(3)==TT(1)) 
                           TrisAux(t,:)=Cell.Tris{c}(t,:);
                           Included(t)=1;
                           break
                    end 
                end 
            end
            if sum(Included)==length(Included)
                TrisAux(1,:)=Cell.Tris{c}(1,:);
                NotAllIncluded=false;
            end
        end
        
        % compute volume     
        auxV=0;
        for t=1:size(TrisAux,1)
            if TrisAux(t,3)<1
                YTri=[Y.DataRow(TrisAux(t,[1 2]),:); Y.DataRow(abs(TrisAux(t,3)),:)];
            else 
                YTri=[Y.DataRow(TrisAux(t,[1 2]),:); Cell.FaceCentres.DataRow(TrisAux(t,3),:)];
            end 
            T=det(YTri)/6;
            auxV=auxV+T;       
        end 
        
        % if the volume is negative switch two the other option
        if auxV<0
            TrisAux=[TrisAux(:,2) TrisAux(:,1) TrisAux(:,3)];
        end 
        
        % Correct Cell and faces Data
        Cell.Tris{c}=TrisAux;
        aux1=1;
        for s=1:Cell.Faces{c}.nFaces
            if Cell.Faces{c}.Vertices{s}(1)==TrisAux(aux1,2) &&...
               Cell.Faces{c}.Vertices{s}(2)==TrisAux(aux1,1)
               Cell.Faces{c}.Vertices{s}=flip(Cell.Faces{c}.Vertices{s});
               if length(Cell.Faces{c}.Vertices{s})==3
                    Cell.Faces{c}.Vertices{s}=abs(TrisAux(aux1,:));
                    Cell.Faces{c}.Tris{s}=TrisAux(aux1,:);
               else 
                  Cell.Faces{c}.Tris{s}=TrisAux(aux1:aux1+length(Cell.Faces{c}.Vertices{s})-1,:);
               end 
               Cell.AllFaces.Vertices{Cell.Faces{c}.FaceCentresID(s)}=Cell.Faces{c}.Vertices{s};
            end 
            if length(Cell.Faces{c}.Vertices{s})==3
                aux1=aux1+1;
            else 
                aux1=aux1+length(Cell.Faces{c}.Vertices{s});
            end 
        end
    end 
    
    % Recompte Volume and Surface Area (This can be improved) 
    if Recompute
        [Cell]=BuildEdges(Cell,Y);
        [Cell]=ComputeCellVolume(Cell,Y);
        Cell.AllFaces=Cell.AllFaces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
        Cell.AllFaces=Cell.AllFaces.ComputeEnergy(Set);
    end 

end


