function createImageJROIsFromLabeledROI(labeledROI,RM)
% Generates a new set of ROIs in ImageJ's ROI Manager from the input
% labeledROI.
%
% Example: createImageJROIsFromLabeledROI(labeledROI,RC)
% 
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

newROIs=[];
uniqueValues=unique(labeledROI(labeledROI>0));
for(i=1:length(uniqueValues))
    boundaryPts = bwboundaries(labeledROI==uniqueValues(i));
    boundaryPts = boundaryPts{1}-1; % subtract one because Java arrays start from 0, rather than 1.
    newROIs = cat(1,newROIs,ij.gui.PolygonRoi(boundaryPts(:,2),boundaryPts(:,1),size(boundaryPts,1),ij.gui.Roi.POLYGON));
end

% Clear existing ROIs
if RM.getCount~=0
    RM.runCommand('Deselect');
    RM.runCommand('Delete');
end

% Add new ROIs
for i=1:length(newROIs)
    RM.addRoi(newROIs(i));
end

end