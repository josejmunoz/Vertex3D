function [] = paintOneCell(idCell)
%PAINTONECELL Summary of this function goes here
%   Detailed explanation goes here

inputDir = '/Users/pablovm/PostDoc/VertexModel3D/Vertex3D/Result/Relevant/NoAblationNoContractility_NoRemodel_S2_0.1_3x3/Analysis/';
matFiles = dir(strcat(inputDir, '*.mat'));

figure

for numFile = 1:size(matFiles, 1)
    load(strcat(inputDir, 'resultingImage_', num2str(numFile),'.mat'))
    resultingImage = imresize3(resultingImage, size(resultingImage).*2, 'nearest');
    [x, y, z] = ind2sub(size(resultingImage), find(resultingImage == idCell));

    shp = alphaShape(x,y,z);
    shp.Alpha = 1;
    plot(shp, 'FaceColor', [57, 156, 203]/255, 'EdgeColor', 'none', 'AmbientStrength', 0.3, 'FaceAlpha', 0.5);
    %hold on;
    %pcshow([x,y,z], [57, 156, 203]/255)
    xlim([0, size(resultingImage, 1)]);
    ylim([0, size(resultingImage, 2)]);
    zlim([0, size(resultingImage, 3)]);
    %axis equal
    camlight left;
    camlight right;
    lighting flat
    material dull
    
    newFig = gca;
    newFig.XGrid = 'off';
    newFig.YGrid = 'off';
    newFig.ZGrid = 'off';
    
    movieInfo(numFile) = getframe();
    %hold off;
end

h = figure,
movie(h),

end

