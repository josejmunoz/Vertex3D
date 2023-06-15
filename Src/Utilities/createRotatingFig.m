function [] = createRotatingFig(fileName)
%CREATEROTATINGFIG Summary of this function goes here
%   Detailed explanation goes here
%
%   You can use it like:
%   files = dir('/Users/pablovm/Dropbox (UCL)/TetsToShowIntercalation/*.fig')
%   for numFile = 1:length(files)
%       uiopen(fullfile(files(numFile).folder, files(numFile).name), 1)
%       createRotatingFig(fullfile(files(numFile).folder, files(numFile).name))
%   end

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

