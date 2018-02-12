function labeledROI = createLabeledROIFromImageJPixels(imgSize,roiObjects)
% Creates a labeled ROI from the ROIs defined in ImageJ's ROI Manager.
%
% Example: labeledROI = createLabeledROIFromImageJPixels([512 512],RC.getRoisAsArray)
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

labeledROI = zeros(imgSize);
nROIs = length(roiObjects);
for(i=1:nROIs)
    % Get center location for ROI object
    X = roiObjects(i).getXBase+1; % add one because MATLAB arrays start at 1, while Java arrays start at 0.
    Y = roiObjects(i).getYBase+1;
    
    % Get local mask for ROI object
    localCellMask = roiObjects(i).getMask();
    height = localCellMask.getHeight();
    width  = localCellMask.getWidth();
    boundedPixels = double(localCellMask.getPixels());
    localCellImg = reshape(boundedPixels,[width,height]);
    localCellImg(localCellImg==-1) = i;
    
    % Add ROI to labeled ROI
    labeledROI(Y+[1:height],X+[1:width]) = localCellImg';
end