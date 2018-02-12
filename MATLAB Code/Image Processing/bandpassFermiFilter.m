function map = bandpassFermiFilter(map,wavelengthLower,wavelengthUpper,spatialResolution,zscoreData,T,ROI)
%function map = bandpassFermiFilter(map,wavelengthLower,wavelengthUpper,spatialResolution,zscoreData,T,ROI)
%
% Performs a bandpass fermi-filter at the desired wavelength cutoff and
% spatial resolution on a complex field image (both in microns). The 
% returned image is the filtered complex field. The complex field map is 
% filtered by performing an individual highpass filter on the respective 
% real and complex values that makeup a complex field.
%
% by David Whitney, Max Planck Florida Institute, 2017
% email: david.whitney@mpfi.org
if(nargin < 5), zscoreData = false;             end
if(nargin < 6), T = 0.05;                       end
if(nargin < 7), ROI = logical(ones(size(map))); end

% Pad Image
imsize = size(map);
[map startingPosition] = padImage(map);
ROI = logical(padImage(ROI));

% Z-Score Image and Remove Any Pixels Outside ROI
if(zscoreData)
    map = map-mean(map(ROI));
    map = map./std(map(ROI));
    map(ROI == 0) = 0;
end

% Bandpass Filter Image
if(wavelengthLower ~= -1) %lowpass filter only
    map = fermifilter(map,wavelengthLower,spatialResolution,T,ROI);
end
if(wavelengthUpper ~= -1) %highpass filter only
    map = map-fermifilter(map,wavelengthUpper,spatialResolution,T,ROI);
end

% Z-Score Image and Remove Any Pixels Outside ROI
if(zscoreData)
    map = map-mean(map(ROI));
    map = map./std(map(ROI));
    map(ROI == 0) = 0;
end

% Crop Image To Original Proportions
map = removePaddingToImage(map,startingPosition,imsize);
ROI = removePaddingToImage(ROI,startingPosition,imsize);

%% Miscellaneous Functions
function im = removePaddingToImage(im,startLocation,imsize)
% function im = cropIm(im,dims)
% takes in the dimensions specified by dims and crops the image around the 
% image's center

im = im(startLocation(1):(startLocation(1)+imsize(1)-1),startLocation(2):(startLocation(2)+imsize(2)-1));

function [paddedIm startingPosition] = padImage(im,paddingValue)
% pads image to nearest value of 2^n and with a padding value defined by paddingValue

if(nargin<2)
    paddingValue = 0;
end

paddedImageSize = round(2.^(ceil(log(max(size(im)))/log(2))));
paddedIm = paddingValue.*ones(paddedImageSize);
startingPosition = round((paddedImageSize-size(im))/2);
for(cdim = 1:2) %Ensure startingPosition is within map
    if(startingPosition(cdim)<1)
        startingPosition(cdim) = 1;
    elseif(startingPosition(cdim)>size(im,cdim))
        startingPosition(cdim) = size(im,cdim);
    end
end
paddedIm(startingPosition(1):(startingPosition(1)+size(im,1)-1),startingPosition(2):(startingPosition(2)+size(im,2)-1)) = im;

function filtim=fermifilter(im,mm_bound,pix_res,T,ROI)
%This is a function in which a spatial fermi filter is applied to an image
%in fourier space.
%
%USAGE
%function filtim=fermifilter(im,mm_boundLower,mm_boundUpper,pix_res,T,ROI)
%
%VARIABLE DEFINITIONS
%im - Input image to be filtered
%mm_bound - The wavelength in µm at which the fermi filter is defined to give a 
%			50% response.  This is effectively the border of the low pass
%			filter to be applied ot the image.
%pix_res - The pixel resolution for the image. 
%filtim - The filtered image to be output.
%T - temperature of filter (should be 0.05)
%ROI - roi for image to be filtered

%From the image size, pixel resolution, and desired 50% filter response, calculate
%the appropriate pixel cutoff for the filter to be applied to the 2D fft. 
imsize=size(im);
ap_size=pix_res*imsize(2);
pix_bound=ap_size./mm_bound;

%Create a matrix that is the same size as the input image and will serve as a map 
%of the distance from the center of the 2D fft image.
[X,Y]=meshgrid(1:imsize(2),1:imsize(1));
X=X-(1+imsize(2)/2);
Y=Y-(1+imsize(1)/2);
dist=sqrt(X.^2+Y.^2);

%Define the fermi filter to be used
fermi2=1./(1+exp(-((pix_bound-abs(dist))./(T*pix_bound))));
fermi2=fftshift(fermi2);

%calculate the 2D fft of the input image and apply the fermi filter to that image.
%next, take the inverse transform in order to get the filtered image back. 
fftedImage  = fft2(im ,imsize(1),imsize(2));
fftedROI    = fft2(ROI,imsize(1),imsize(2));
filtim      = ifft2(fftedImage.*fermi2);
filtROI     = ifft2(  fftedROI.*fermi2);
filtim(ROI) = filtim(ROI) ./ filtROI(ROI);