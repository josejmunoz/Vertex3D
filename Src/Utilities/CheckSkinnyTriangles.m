function [s]=CheckSkinnyTriangles(Y1,Y2,cellCentre)
	YY12=norm(Y1-Y2);
	Y1=norm(Y1-cellCentre);
	Y2=norm(Y2-cellCentre);
	
	if YY12>2*Y1 || YY12>Y2*2
    	s=true;
	else
    	s=false;
	end
end 