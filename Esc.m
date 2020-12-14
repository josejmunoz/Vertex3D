function esc=Esc()
% Returns separator for folders
machine=computer;
if strcmp(machine(1:3),'PCW')==1
    esc='\';
% Cluster or MAC    
elseif strcmp(machine(1:3),'MAC')==1 || strcmp(machine(1:3),'GLN')==1
    esc='/';
else
    esc='/';
end
end