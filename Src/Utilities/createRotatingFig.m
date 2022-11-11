function [] = createRotatingFig(fileName)
%CREATEROTATINGFIG Summary of this function goes here
%   Detailed explanation goes here
v = VideoWriter(fileName, 'MPEG-4'); %% if linux use other than mp4
v.Quality = 30;
open(v);

axis equal
%lighting flat
material dull

axis off

for k = 1:360
    camorbit(1,0,'data',[0 0 1])
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
end

