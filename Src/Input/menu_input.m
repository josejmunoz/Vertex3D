function [Geo,Set] = menu_input(inputMode, batchMode)
%MENU_INPUT Summary of this function goes here
%   Detailed explanation goes here
switch inputMode
    case 1
        Stretch
        if ~batchMode
            disp('STRECH SIMULATION');
        end
    case 2
        StretchBulk
        if ~batchMode
            disp('STRECH BULK SIMULATION');
        end
    case 3
        Compress
        if ~batchMode
            disp('COMPRESSION SIMULATION');
        end
    case 4
        Remodelling_Bubbles
        if ~batchMode
            disp('REMODELLING WITH BUBBLES SIMULATION');
        end
    case 5
        Remodelling_Voronoi
        if ~batchMode
            disp('REMODELLING WITH VORONOI SIMULATION');
        end
    case 6
        NoBulk
        if ~batchMode
            disp('REMODELLING WITHOUT BULK');
        end
    case 7
        NoBulk_110
        if ~batchMode
            disp('REMODELLING WITHOUT BULK 110 CELLS');
        end
    otherwise
        error('Incorrect mode selected');
end
end

