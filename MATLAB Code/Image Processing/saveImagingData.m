function saveImagingData(tifStack,saveDirectory,useImageRegistration,chunkSize)
% saves tif stack to the specified saveDirectory via ImageJ.
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

if(nargin<3), useImageRegistration = true; end % specifies whether the saved images were registered. used to denote the tag used for saving
if(nargin<4), chunkSize = 1000; end % number of imaging frames saved in each multipage tif

% Save registered stack
numberOfFrames = size(tifStack,3);
numberOfChunks = ceil(numberOfFrames/chunkSize);
for(ii=1:numberOfChunks)
    % Send stack to MIJ
    selectedFrames = [1:chunkSize]+(ii-1)*chunkSize;
    selectedFrames = selectedFrames(selectedFrames<numberOfFrames);
    MIJImgStack = MIJ.createImage('ReferenceStack',single(tifStack(:,:,selectedFrames)),true);

    % Save tif stack in MIJ
    if(useImageRegistration), tag = 'registered'; else, tag = 'unregistered'; end
    fileName = strrep(sprintf('%s%s_%05d.tif',saveDirectory,tag,ii),'\','\\');
    ij.IJ.run(MIJImgStack,'Save', sprintf('Tiff, path=[%s]',fileName));
    MIJImgStack.close();
end