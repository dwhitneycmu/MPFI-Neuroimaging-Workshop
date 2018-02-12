function tifStack = readMultipageTifFiles(folderLocation)
% reads multipage tif stacks from the specified folder location into a
% single 3d imaging stack.
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

if(nargin<1), folderLocation  = 'E:\WinterCourse\RAW Data\RAW\'; end % Defines where our raw tif stacks are saved.

tifStack = [];
tifFiles = dir([folderLocation '*.tif']);
for(currentFile = 1:length(tifFiles))
    disp(['Reading Imaging Stack ' num2str(currentFile) ' Out Of ' num2str(length(tifFiles))]);

    % Specify stack name
    filePath = [folderLocation tifFiles(currentFile).name];

    % Read images into tifStack
    tifStack = cat(3,tifStack,read_Tiffs(filePath,1,50));
end