function imgStack = shiftImageStack(imgStack,regOffsetsYX)
% imgStack = shiftImageStack(imgStack,regOffsetsYX)
% Shift each image in the input image stack using the specified (Y,X) registration offsets
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

for(d = 1:size(imgStack,3))   
    try
        imgStack(:,:,d) = imtranslate(imgStack(:,:,d),regOffsetsYX(d,:));
    catch
        disp('Error')
    end
end