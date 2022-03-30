function CellStr = BuildStructArray(n, fields)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BuildStructArray:										  
	%   Creates an array of n structs, each with fields. MATLAB imposes 
	%   that all structs have the same fields on creation. That is, if 
	%   we were to append an empty struct to CellStr, MATLAB errors out.
	% Input:															  
	%   n       : Length of the struct array  
	%   fields  : Fields each struct of the array will have
	% Output:															  
	%   CellStr : Array of structs
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	CellStr = struct();
	for f = 1:length(fields)
		CellStr.(fields(f)) = {};
	end
	for c = 2:n
		temp_str = struct();
		for f = 1:length(fields)
			temp_str.(fields(f)) = {};
		end
		CellStr(c) = temp_str;
	end
end