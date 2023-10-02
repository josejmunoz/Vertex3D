function v=ComputeCellVolume(Cell)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ComputeCellVolume:										  
	%   Computes Cell Volume
	% Input:															  
	%   Cell : Cell object for which the volume is calculated			  
	% Output:															  
	%   v : Cell volume    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	v = 0;
    for f = 1:length(Cell.Faces)
		face = Cell.Faces(f);
		for t=1:length(face.Tris)
			y1 = Cell.Y(face.Tris(t).Edge(1),:) - Cell.X;
			y2 = Cell.Y(face.Tris(t).Edge(2),:) - Cell.X;
            y3 = face.Centre - Cell.X;
			Ytri = [y1; y2; y3];
            
            currentV = det(Ytri)/6;
            % if the volume is negative switch two the other option
            if currentV < 0
                Ytri = [y2; y1; y3];
                currentV = det(Ytri)/6;
            end

			v = v + currentV;
		end
    end

%     if v < 0
%         v = 0;
%         for f = 1:length(Cell.Faces)
%             face = Cell.Faces(f);
%             for t=1:length(face.Tris)
%                 y1 = Cell.Y(face.Tris(t).Edge(1),:);
%                 y2 = Cell.Y(face.Tris(t).Edge(2),:);
%                 y3 = face.Centre;
%                 Ytri = [y2; y1; y3];
%                 currentV = det(Ytri)/6;
% %                 % if the volume is negative switch two the other option
% %                 if currentV < 0
% %                     Ytri = [y2; y1; y3];
% %                     currentV = det(Ytri)/6;
% %                 end
%                 
%                 v = v + currentV;
%             end
%         end
%     end
end 