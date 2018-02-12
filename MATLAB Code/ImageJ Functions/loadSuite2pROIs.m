% Loads ROIs generated from Suite2p, and sends them to ImageJ
load([baseDirectory 'analyzedData\Suite2p\F_FerretImagingCourse_2016-01-14_plane1.mat'],'ops','stat')

% Generated a labeled ROI of the Suite2pROIs
ROI      = zeros(ops.Ly,ops.Lx);
validROI = zeros(length(ops.yrange),length(ops.xrange));
nCells = length(stat);
LUT=hsv(nCells); rng(1);
LUT=LUT(randperm(nCells),:);
for(i=1:nCells)
    validROI(stat(i).ipix)=i;
end
ROI(ops.yrange(1)-1+[1:length(ops.yrange)],...
    ops.xrange(1)-1+[1:length(ops.xrange)])=validROI;

% Optionally display in MATLAB the Suite2p ROIs
showSuite2pROIs = false;
if(showSuite2pROIs)
    ROI_RGB = repmat(ROI>0,[1 1 3]).*ind2rgb(ROI,LUT);
    figure; imagesc(ROI_RGB); colormap(hsv); axis image; axis off;
end

% Send ROIs to ImageJ
createImageJROIsFromLabeledROI(ROI,RC)