function tifStack = read_Tiffs(filePath,imgScaling,updateFrequency,dataType,useWaitBar,codeVersion)
% tifStack = read_Tiffs(filePath,updateFrequency)
%
% A fast method to read a tif stack: calls the tiff library directly.
% Used the same method described in this web-link: 
% http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/
% 
% filePath - fileDirectory and name of image stack
% imgScaling - Scales size of images by specified scaling factor (default is no scaling).
%              0.5 drops image size by half, while 2x doubles image size
% updateFrequency - Based on the number of images, how often in percentage
% of files read to inform user of current progress (0-100%)
% dataType - data type of the imaging data (default is 'uint16', but can be 'single' or 'float')
% useWaitBar - Displays the progress of reading in the current stack.
% codeVersion - Either uses the tifflib function to read images (older method), or
%               uses the newer Tiff MATLAB object to read images (new, default method).
%
% Written by David Whitney (10/3/2013)
% Updated on 8/1/2017
% David.Whitney@mpfi.org
% Max Planck Florida Institude

if(nargin<2), imgScaling      = 0.5;      end
if(nargin<3), updateFrequency = 100;      end
if(nargin<4), dataType        = 'uint16'; end
if(nargin<5), useWaitBar      = false;    end
if(nargin<6), codeVersion     = 'Tiff';   end
tic;

disp(['Reading Image Stack - ' filePath]);

% Read TIF Header and Setup Information
InfoImage=imfinfo(filePath);
    xImage=InfoImage(1).Width;
    yImage=InfoImage(1).Height;
    NumberOfImages=length(InfoImage);
    disp(['Finished Reading Image Header - ' num2str(toc) ' seconds Elapsed']);

% setup wait bar
if(useWaitBar)
    h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
    currentPosition = get(h,'Position');
    offset = 100;
    set(h,'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
    currentPosition = get(get(h,'Children'),'Position');
    set(get(get(h,'Children')),'Position',[currentPosition(1)-offset/2 currentPosition(2) currentPosition(3)+offset currentPosition(4)]);
else
    updateFrequency = round((updateFrequency/100)*NumberOfImages);
end

% Initialize MATLAB array to contain tif stack
scaledX = round(xImage*imgScaling);
scaledY = round(yImage*imgScaling);
tifStack     = zeros(scaledY,scaledX,NumberOfImages,dataType);
currentImage = zeros(yImage,xImage);
   
% Get imaging data
switch codeVersion
    case 'tifflib'
        % uses the tifflib function to read images fast
        FileID = tifflib('open',filePath,'r');
        rps    = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
        for i=1:(NumberOfImages)
            if(useWaitBar)
                waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
            else
                if(mod(i+1,updateFrequency)==0)
                    disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                end
            end

            tifflib('setDirectory',FileID,i-1);
            warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
            % Go through each strip of data.
            rps = min(rps,yImage);
            for r = 1:rps:yImage
              row_inds = r:min(yImage,r+rps-1);
              stripNum = tifflib('computeStrip',FileID,r)-1;
              currentImage(row_inds,:) = tifflib('readEncodedStrip',FileID,stripNum);
            end
            if(imgScaling ~= 1 && imgScaling>0)
                tifStack(:,:,i) = imresize(currentImage,[scaledY scaledX]); % Scales image size
            else
                tifStack(:,:,i) = currentImage;
            end
        end
        tifflib('close',FileID);
        disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
        warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
    case 'Tiff'
        % Setup TIF object and Read-In Basic Information
        hTif = Tiff(filePath);
          
        warning('OFF','MATLAB:imagesci:tiffmexutils:libtiffWarning');
        for i=1:(NumberOfImages)
            if(useWaitBar)
                waitbar(i/NumberOfImages,h, ['Image ' num2str(i) ' of ' num2str(NumberOfImages) ' - ' num2str(toc) 's Elapsed - ' num2str((NumberOfImages-i)*toc/i) 's Left']);
            else
                if(mod(i+1,updateFrequency)==0)
                    disp([num2str(round(100*i/NumberOfImages)) '% Done Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
                end
            end

            if(imgScaling ~= 1 && imgScaling>0)
                tifStack(:,:,i) = imresize(hTif.read(),[scaledY scaledX]); % Scales image size
            else
                tifStack(:,:,i) = hTif.read();
            end
            if(i == NumberOfImages)
                hTif.close();
                warning('ON','MATLAB:imagesci:tiffmexutils:libtiffWarning')
                disp(['Finished Reading Image Stack - ' num2str(toc) ' seconds Elapsed']);
            else
                hTif.nextDirectory();                
            end
        end
end

if(useWaitBar),close(h); end