function neuropilRois = generateNeuropilROIs(rois,diskRadius,dilationRadius)
% takes cellular imageJ ROIs and creates neuropil masks centered around each roi.
% The neuropil masks are disks projected away from each cell with a radius equal
% to the pixel size specified by diskRadius. Neuropil masks will not
% include any pixels contained in any of the cellular rois. 
%
% rois - can either be a path to the location of ROIs or are the cellular rois
%        passed into MATLAB with rm.getRoisAsArray
% diskRadius - size of the neuropil masks
% dilationRadius - number of pixels to dilate cellular masks by (0 is default)
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.
if(nargin<2), diskRadius = 15;    end
if(nargin<3), dilationRadius = 1; end

% Setup ROI Manager with MIJ
RC = ij.plugin.frame.RoiManager();
RM = RC.getInstance();

% Get cellular ROIs
if(ischar(rois)) % if a path is specified we'll load ROIs
    % Specify ROI path
    ROIPath = sprintf(rois);

    % Get cellular ROIs as arrays
    RM.runCommand('Open',ROIPath)
    rois = RM.getRoisAsArray;
end

%% Step 0: Pre-Process ROIs and explicitly convert them to a ShapeRoi. 
% Sometimes ImageJ will use a different ROI type, and this code will throw
% an error as it expects the ROIs to inherit specific methods from the 
% ShapeROI class.
cellRois = rois;
for i =1:length(cellRois)
    cellName = cellRois(i).getName;
    cellRois(i) = ij.gui.ShapeRoi(cellRois(i));
    cellRois(i).setName(cellName);
end

%% Step 1: Dilate cellular ROIs
emptyOverlay=ij.gui.Overlay; % allows us to put an empty here
dilateCellularROIs=emptyOverlay.toArray;
for i=1:length(cellRois)
    dilateCellularROIs(i) = ij.plugin.RoiEnlarger.enlarge(cellRois(i).clone,dilationRadius); % dilates cellular roi
    dilateCellularROIs(i) = ij.gui.ShapeRoi(dilateCellularROIs(i)); % convers from a traced ROI to a composite roi (needed to get logical operations like OR or XOR)
    dilateCellularROIs(i).setName(cellRois(i).getName);
    %dilateCellularROIs(i).setName([cellRois(i).getName.toCharArray' '_Test']); %only used for debugging purposes
end

%% Create a mask set where we project discs away from the centers of each cellular ROI
% Step 2: Create a ROI that contains regions located within our cell ROIs. All future 
% neuropil masks will exclude this region.
excludedRegions=dilateCellularROIs(1).clone;
for i=2:length(cellRois)
    excludedRegions=excludedRegions.or(dilateCellularROIs(i));
end

% Step 3: Generate a neuropil mask by creating a disk roi centered around each cell. However, 
% we will exclude all regions that include our cellular ROIs.
neuropilRois = emptyOverlay.toArray;
for i=1:length(cellRois)
    % Create new oval roi
    cellPosition = cellRois(i).getBoundingRect();
    x = cellPosition.getX+0.5*cellPosition.getWidth;
    y = cellPosition.getY+0.5*cellPosition.getHeight;
    diskDiameter = 2*diskRadius+1;
    ovalROI  = ij.gui.OvalRoi(x-diskRadius,y-diskRadius,diskDiameter,diskDiameter); % start with an oval roi
    shapeROI = ij.gui.ShapeRoi(ovalROI); % convert oval roi to a shape roi    
    
    % Save the new roi, but we will exclude pixels contained within
    % cellular ROIs
    neuropilRois(i) = shapeROI.not(excludedRegions);
    
    % Rename the neuropil ROIs
    cellName = char(cellRois(i).getName);
    neuropilRois(i).setName(strrep(cellName,'Cell','Neuropil'));
end    

% Step 4: Load neuropil masks into imageJ
for i=1:length(neuropilRois)
    RM.addRoi(neuropilRois(i));
end
return