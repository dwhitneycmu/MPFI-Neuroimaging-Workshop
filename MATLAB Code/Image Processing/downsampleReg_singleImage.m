function [regOffsetsXY,varargout]=downsampleReg_singleImage(sourceImg,template,downsampleRates,maxMovement)
%  [regOffsetsXY imout]=downsampleReg_singleImage(imgStack)
% 
% DownsampleReg: 
% Original implementation by Theo Walker (Jan. 3, 2014)
% Revised by David Whitney (Aug. 4, 2017)

% === Description === %

%This registers time-series data, removing motion artifacts. 
%It is assumed that your input directory
%contains a series of TIF images; an output directory of registered TIFs is
%produced.

% === Algorithm Explanation === %

%So, the naive way to do register a pair of images is to consider all 
%possible shifts of one image relative to the other, and find where they 
%correlate best. This algorithm is an improvement on that idea.

%The way this works is super simple. What we're going to do is downsample
%the image to something _very_ small, say 1/16 of the original size, and
%register all the slices using correlation.
%Then based on the 1/16 registration, we do a finer registration (1/8) that is
%constrained by the previous registration result, and so on. The constraining
%helps keep the algorithm from focusing on small details that don't matter 
%(i.e., noise) and also speeds up the runtime.

% === Parameters === %

%Downsample rates: As described in algorithm. Examples follow.
%
%downsampleRates = [1/16, 1/8, 1/4, 1/2, 1] is the default. 
%
%downsampleRates = [1/16, 1/4, 1] runs slightly faster. It does not
%sacrifice any theoretical accuracy, but it will be confused by noise a bit
%easier than the default will.
%
%downsampleRates = [1/16] will do a registration _just_ based on the 16x
%downsampled image. It will run pretty quick, but the output will only be
%accurate to +/- 8 pixels; it won't do fine motion correction.

%Max movement: The greatest degree of motion correction allowed. 
%e.g, if this is at 1/4 and your image is 512x512, this will fix movements 
%of up to 128 pixels.
%Makes runtime of the first iteration a little shorter to set this lower.

%  === Code begins here! === %


%Thresholding is a nice preprocessing step, improves the results a bit.
%imgStack(find(imgStack < mean(mean(mean(imgStack))))) = mean(mean(mean(imgStack)));

% Setup constants
if(nargin<3), downsampleRates = [1/16, 1/8, 1/4, 1/2, 1]; end
if(nargin<4), maxMovement     = 1/8;                      end %1/4;
[height,width] = size(sourceImg);
regOffsetsYX   = zeros(2,1); %stores the output of the registration (Y,X format)

% Loop through each downsample step until we get a progressively better
% estimate of the X,Y shift
for r=1:length(downsampleRates)
    % downsample the images. Bilinear interpolation's going to be good here. 
    % Default is bicubic, which is both slow and overkill.
    %disp(['registering images, iteration ' num2str(r)]);
    downHeight  = height*downsampleRates(r);
    downWidth   = width*downsampleRates(r);
    regImg      = imresize(sourceImg,[downHeight,downWidth],'bilinear');
    templateImg = imresize(template, [downHeight,downWidth],'bilinear');
        
    % determine the X and Y shifts needed to register the source image to the template image;
    if r==1
        %we're finding the initial offset
        minOffsetY = -round(maxMovement*downHeight/2);
        maxOffsetY =  round(maxMovement*downHeight/2);
        minOffsetX = -round(maxMovement*downWidth/2);
        maxOffsetX =  round(maxMovement*downWidth/2);

    else
        %we are refining an earlier offset
        minOffsetY = regOffsetsYX(1)*downsampleRates(r) - downsampleRates(r)/downsampleRates(r-1)/2;
        maxOffsetY = regOffsetsYX(1)*downsampleRates(r) + downsampleRates(r)/downsampleRates(r-1)/2;
        minOffsetX = regOffsetsYX(2)*downsampleRates(r) - downsampleRates(r)/downsampleRates(r-1)/2;
        maxOffsetX = regOffsetsYX(2)*downsampleRates(r) + downsampleRates(r)/downsampleRates(r-1)/2;
    end
    bestCorrValue = -1;
    bestCorrX = 0;
    bestCorrY = 0;
    for y=minOffsetY:maxOffsetY
        for x=minOffsetX:maxOffsetX
            %determine the offsets in X and Y for which the overlap between the images correlates best
            subTemplateY1 = 1+max(y,0);
            subTemplateY2 = downHeight+min(y,0);
            subTemplateX1 = 1+max(x,0);
            subTemplateX2 = downWidth+min(x,0);

            subRegY1 = 1+max(-y,0);
            subRegY2 = downHeight+min(-y,0);
            subRegX1 = 1+max(-x,0);
            subRegX2 = downWidth+min(-x,0);

            subTemplateImg = templateImg(subTemplateY1:subTemplateY2,subTemplateX1:subTemplateX2);
            subRegImg = regImg(subRegY1:subRegY2,subRegX1:subRegX2);

            corrValue = corr2(subRegImg, subTemplateImg);
            if corrValue > bestCorrValue
                bestCorrX = x;
                bestCorrY = y;
                bestCorrValue = corrValue;
            end
        end
    end
    regOffsetsYX = [bestCorrY bestCorrX]*1/downsampleRates(r);
end

% now we have the regOffsets; apply them to the images and output results.
if nargout>1
    varargout{1}=shiftImageStack(sourceImg,regOffsetsYX);
end
regOffsetsXY=(regOffsetsYX(:,[2 1])); % shift into x,y format