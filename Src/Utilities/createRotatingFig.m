function [] = createRotatingFig(fileName)
%CREATEROTATINGFIG Summary of this function goes here
%   Detailed explanation goes here
v = VideoWriter(strcat(fileName, '.mp4'), 'MPEG-4');
open(v);
v.Quality = 40;

axis equal
camlight left;
camlight right;
lightning flat
material dull

newFig = gca;
newFig.XGrid = 'off';
newFig.YGrid = 'off';
newFig.ZGrid = 'off';

for k = 1:360
    camorbit(1,0,'data',[0 0 1])
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
end

