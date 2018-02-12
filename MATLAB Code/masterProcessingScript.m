%% Winter 2018 Imaging Course Script (by David Whitney)
%% Script processing steps
loadImagingData = true;  % Loads imaging data into MATLAB
    registerImagingStack = false;  % This flag denotes whether we want to register our full imaging stack. Else it loads previously registered data.
    useDownsampledStack  = true;   % This flag can be used if you want just to use the 200X temporally downsampled stack. registerImagingStack flag must be set to false.
    ROIsToLoad = 'CellMagicWand';  % ROIs can be either 'None','CellMagicWand','Suite2p'
recomputeTraces      = false;   % Recomputes cellular responses (can take awhile)
    neuropilCorrectionType = 'none'; % can be 'none', fixed' (i.e. 0.7), and 'adaptive'
showCellsVsNeuropil  = false;   % Compares the correlation of cells and neuropil responses
showDeltaFOverF      = false;   % Shows how the percentile filter works
showTuningCurves     = false;   % Shows tuning curves and trial evoked responses
showFunctionalMaps   = false;   % Shows functional maps

%% Specify where you saved your course data and change MATLAB course directory. We also need to specify where MATLAB and ImageJ are
baseDirectory = 'H:\WinterCourse\';
ImageJ_Path   = 'H:\WinterCourse\ImageJ\Fiji.app\';
MATLAB_Path   = 'C:\Program Files\MATLAB\R2017b\';

%% Add functions to path, connect MATLAB with MIJ, and then setup ROI Manager
addpath(genpath([baseDirectory '\MATLAB']))
addpath([ImageJ_Path '\scripts\']);
Miji
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();

%% Load TIF stacks into MIJ and select ROIs with Cell Magic Wand tool
% Open up ImageJ, load your image stack, and select or load your ROIs
if loadImagingData    
    % Loads registered imaging sequence into MIJ
    tifFilesLocation = [baseDirectory 'RAW Data\RAW\'];        % Defines where our raw tif stacks are located
    saveLocation     = [baseDirectory 'RAW Data\Registered\']; % Defines where our registered imaging stacks are saved
    spatialResolution = 1.1375; % Microns per pixel
    if(registerImagingStack)    
        % Generate imaging stack from multipage tif files
        tifStack = readMultipageTifFiles(tifFilesLocation);

        % Register imaging data
        [tifStack,xyShifts] = imageRegistration(tifStack);

        % Save registered imaging data and trigers
        saveImagingData(tifStack,saveLocation,true);
    else
        if(useDownsampledStack)
            tifStack = read_Tiffs([baseDirectory 'RAW Data\rawData_200xDownsampled.tif'],1,100);
        else
            tifStack = readMultipageTifFiles(saveLocation);
        end
        MIJImgStack = MIJ.createImage('Imaging data',tifStack,true); 
    end
    
    % Loads cellular ROIs into the ROI Manager
    if RC.getCount~=0 %Clears any ROIs already in memory
        RC.runCommand('Deselect');
        RC.runCommand('Delete');
    end
    switch ROIsToLoad
        case 'None'
            % Don't load any ROIs
        case 'CellMagicWand'
            RC.runCommand('Open', strrep([baseDirectory 'analyzedData\RoiSet.zip'],'\','\\')); 
        case 'Suite2p'
            loadSuite2pROIs;
    end
        
    % Generate neuropil ROIs
    cellNumber = RC.getCount();     % You should have 199 if you use the ROIs provided
    %generateNeuropilROIs(RC.getRoisAsArray,round(30/spatialResolution),round(5/spatialResolution)); % Generates neuropil ROIs - Adds an an additional 199 ROIs
    generateNeuropilROIs(RC.getRoisAsArray); % Generates neuropil ROIs - Adds an an additional 199 ROIs
end

%% Extract fluorsescence timecourse for each cellular ROI using MIJ
cd([baseDirectory 'MATLAB\']);  % Make sure to run this
currentImgStack = ij.IJ.getImage();
if recomputeTraces
    cells = [];
    cells.rawF = [];
    cells.rawF_neuropil = [];
    cells.xPos = zeros(cellNumber,1);
    cells.yPos = zeros(cellNumber,1);
    for i = 1:cellNumber
        % Select cell ROI in ImageJ
        fprintf('Processing Cell %d\n',i)

        % Get cell ROI name and parse out (X,Y) coordinates
        RC.select(i-1); % Select current cell
        [tempLoc1,tempLoc2] = strtok(char(RC.getName(i-1)),'-');
        cells.xPos(i) =  str2double(tempLoc1);
        cells.yPos(i) = -str2double(tempLoc2);

        % Get the fluorescence timecourse for the cell and neuropol ROI by
        % using ImageJ's "z-axis profile" function. We can also rename the
        % ROIs for easier identification.
        for isNeuropilROI = 0:1
            ij.IJ.getInstance().toFront();
            MIJ.run('Plot Z-axis Profile'); % For each image, this outputs four summary metrics: number of pixels (in roi), mean ROI value, min ROI value, and max ROI value
            RT = MIJ.getResultsTable();
            MIJ.run('Clear Results');
            MIJ.run('Close','');
            if isNeuropilROI
                %RC.setName(sprintf('Neuropil ROI %d',i));
                cells.rawF_neuropil(i,:) = RT(:,2);
            else
                %RC.setName(sprintf('Cell ROI %d',i));
                cells.rawF(i,:) = RT(:,2);
                RC.select((i-1)+cellNumber); % Now select the associated neuropil ROI
            end
        end
    end
    
    % Now load our stimulus triggers 
    [cells.stimOnsetTime,cells.stimID,cells.frameTimes]=getStimulusTriggers([baseDirectory 'RAW Data\triggers\']);
    if(useDownsampledStack)
        downsamplingFactor  = 200;
    else
        downsamplingFactor  = 2; 
    end
    cells.frameTimes = cells.frameTimes(1:downsamplingFactor:end); % SPECIAL FLAG - Only needed here because we temporally downsampled the stack by 10x
    
    % Compute the neuropil-contributed signal from our cells, and eliminate
    % them from our raw cellular trace
    for(i = 1:cellNumber)
        % Determine scaling factor
        switch neuropilCorrectionType
            case 'adaptive'
                X = cat(1,cells.rawF_neuropil(i,:),cells.rawF(i,:));
                coeffs=robustfit(X(1,:),X(2,:));
                coeffs(coeffs>1)=1;
                coeffs(coeffs<-1)=-1;
                cells.subtractionFactor(i) = coeffs(2);
            case 'fixed'
                cells.subtractionFactor(i) = 0.7; %0.3-0.6--> Kerlin et al., 2010 (Reid),  0.7--> Chen et al., 2013 (Svoboda), 0.5-0.8--> Golstein et al., 2013 (Pennartz), % 1.0-->Miller et al., 2014 (Yuste lab)
            case 'none'
                cells.subtractionFactor(i) = 0;
        end
        
        cells.rawF(i,:) = cells.rawF(i,:)-cells.subtractionFactor(i)*cells.rawF_neuropil(i,:);
    end

    %% Define a moving and fixed fluorescence baseline using a percentile filter
    cells.rate      = 1/median(diff(cells.frameTimes)); % frames per second
    cells.baseline  = 0*cells.rawF;
    cells.f0        = zeros(cellNumber,1);
    for i =1:cellNumber
        disp(sprintf('Computing baseline for cell %d', i))

        % Compute a moving baseline with a 60s percentile lowpass filter smoothed by a 60s Butterworth filter 
        percentileFiltCutOff = 10;
        lowPassFiltCutOff    = 60; %in seconds
        cells.baseline(i,:)  = baselinePercentileFilter(cells.rawF(i,:)',cells.rate,lowPassFiltCutOff,percentileFiltCutOff);

        % Compute a fixed baseline
        percentileCutOff = 5;
        cells.f0(i) = prctile(cells.rawF(i,:),percentileCutOff,2);
    end
else
    load([baseDirectory 'analyzedData\cellularData_Downsampling2x_NoNeuropilSubtraction.mat']);
end

% Show the response timecourse of all cells and surrounding neuropil
if(showCellsVsNeuropil)
    figure; 
    clippingValue = prctile([cells.rawF(:); cells.rawF_neuropil(:)],99);
    featureType = [1,2,4]; % 1 - cells, 2 - neuropil, 3 - difference
    for ii = 1:length(featureType)
        switch featureType(ii)
            case 1
                feature = 'Cellular';
                img     = cells.rawF;
            case 2
                feature = 'Neuropil';
                img     = cells.rawF_neuropil;
            case 3
                feature = 'Difference';
                img     = cells.rawF - cells.rawF_neuropil;
            case 4
                feature = 'Combined';
                img     = cat(1,cells.rawF,cells.rawF_neuropil);
        end

        % Compute mean correlation between cells or neuropil
        traces = zscore(img,0,2);
        correlationTable = (traces*traces')/size(traces,2);
        meanCorrelation  =  mean(correlationTable(triu(correlationTable,1)~=0));
        if(featureType(ii)==4)
            corrCellsWithNeuropil = mean(diag(correlationTable(1:cellNumber,(cellNumber+1):end)))
        end

        % Show response timecourses as a heat map
        subplot(length(featureType),1,ii); 
            imagesc(img); 
            title([feature ' timecourse: r=' num2str(meanCorrelation)]); 
            ylabel('Raw fluorescence'); 
            xlabel('Time (in frames)');
            caxis([0 clippingValue]);
    end
end
 
%% Compare the moving and fixed baseline measurements
if(showDeltaFOverF)
    figure;
    for i=1:10:cellNumber % Looking at every 10th cell
        hold on; 
        plot(cells.rawF(i,:),'k');
        plot(cells.baseline(i,:),'r');
        plot(cells.f0(i)+0*cells.rawF(i,:),'b');
        legend({'Raw cell trace','Moving filter baseline','Fixed baseline'});
        title(['Cell number ', num2str(i)]);
        ylabel('Raw fluorescence'); 
        xlabel('Time (in frames)');
        pause;
        clf;
    end
end
 
%% Convert the fluorescence timecourses to DeltaF/F using the moving baseline
cells.df = (cells.rawF-cells.baseline)./cells.baseline;
if(showDeltaFOverF)
    figure; imagesc(cells.df);
    title('\DeltaF/F fluorescence timecourse')
    ylabel('\DeltaF/F');
    xlabel('Time (in frames)');
    colorbar;
    caxis([0 3]);
end

% Find stimulus onsets
cells.stimDuration   = 3; % Duration that stimulus was presented (in seconds)
cells.preTrialTime   = cells.stimDuration/2; % Duration of the prestimulus interval (in seconds)
cells.postTrialTime  = cells.stimDuration;   % Duration of the prestimulus interval (in seconds)
cells.uniqStims      = length(unique(cells.stimID));
cells.numberOfStims  = length(cells.stimOnsetTime);
cells.stimStartIndex = zeros(cellNumber,1);
cells.stimStopIndex  = zeros(cellNumber,1);
for i = 1:cells.numberOfStims
    cells.stimStartIndex(i) = find(cells.frameTimes>cells.stimOnsetTime(i),1);
    cells.stimStopIndex(i)  = floor(cells.stimStartIndex(i)+cells.stimDuration*cells.rate);
end
 
% Determine trial-evoked responses 
cells.stimResponse     = [];
cells.meanStimResponse = [];
for stimulus = 1:cells.uniqStims
    stimIndices = find(cells.stimID==stimulus);
    for trialNumber =  1:length(stimIndices)
        for i = 1:cellNumber
            % Determine frames when stimulus was presented
            stimIndex  = stimIndices(trialNumber);
            offsetPre  = round(cells.rate*cells.preTrialTime);
            offsetPost = round(cells.rate*cells.postTrialTime);
            selectedFramesTrace = (cells.stimStartIndex(stimIndex)-offsetPre):(cells.stimStopIndex(stimIndex)+offsetPost); % Here we are padding the response timecourse with some extra pre-trial and post-trial frames to see the rise/decay of the calcium
            selectedFramesMean  = cells.stimStartIndex(stimIndex):cells.stimStopIndex(stimIndex); % Here we only consider responses when the stimulus was presented 
            
            % Extract the trial-evoked timecourse and mean trial response
            cells.stimResponse(i,stimulus,trialNumber,:)   = cells.df(i,selectedFramesTrace);
            cells.meanStimResponse(i,stimulus,trialNumber) = mean(cells.df(i,selectedFramesMean));
        end
    end
end
 
%% Compare different visually evoked responses
if(showTuningCurves)
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = [2 38 69 86] %1:cellNumber
        % Compute summary stats for responses
        n     = size(cells.meanStimResponse,3);
        x     = linspace(0,360,cells.uniqStims);
        y     = squeeze(cells.meanStimResponse(i,:,:));
        yMean = mean(y,2);
        yStd  = std(y,[],2);
        ySEM  = yStd/sqrt(n);
        preferredStimulus = find(yMean(1:(end-1)) == max(yMean(1:(end-1))));
        responseThreshold = yMean(end)+2*yStd(end); % The average response at the preferred stimulus must be 2 standard deviations above the blank condition mean
        IsRespSignificant = yMean(preferredStimulus)>responseThreshold;

        % Compute stats for preferred response
        t = ([1:size(cells.stimResponse,4)]/cells.rate)-cells.preTrialTime;
        yResponse     = squeeze(cells.stimResponse(i,preferredStimulus,:,:));
        yResponseMean = mean(yResponse,1);
        yResponseSEM  = std(yResponse,[],1)/sqrt(size(yResponse,1));

        % Compute stats for blank stimulus (last ID)
        yBlank     = squeeze(cells.stimResponse(i,end ,:,:));
        yBlankMean = mean(yBlank,1);
        yBlankSEM  = std(yBlank,[],1)/sqrt(size(yBlank,1));

        % Show response timcourse for preferred response
        subplot(1,3,1); 
            plot(t,yResponseMean,'-r','lineWidth',3);
            hold on;
            plot(t,yResponse,'--k','Color',0.25*[1,0,0]);
            legend({'Average response','Trial responses'},'Location','northwest');
            xlim([min(t) max(t)]);
            set(gca,'Box','off','XTick',cells.stimDuration*[0 1 2]);
            title(sprintf('Preferred response at %d%s for cell %d',round(360*((preferredStimulus-1)/(length(yMean)-1))),char(176),i));
            ylabel('\DeltaF/F')
            xlabel('Time (seconds)')
            axis square;

        % Show response timecourse with error bars for preferred and blank response
        subplot(1,3,2); 
            dataToShow = 1:1:length(t);
            errorbar(t(dataToShow),yResponseMean(dataToShow),yResponseSEM(dataToShow),'-r','lineWidth',1);
            hold on;
            errorbar(t(dataToShow),yBlankMean(dataToShow),yBlankSEM(dataToShow),'-k','lineWidth',1);
            plot(t(dataToShow),responseThreshold+0*t(dataToShow),'--k');
            legend({'Preferred direction','Blank','Response threshold'},'Location','northwest');
            xlim([min(t) max(t)]);
            set(gca,'Box','off','XTick',cells.stimDuration*[0 1 2]);
            if IsRespSignificant
                title('Preferred response is significant');
            else
                title('Preferred response not significant');
            end
            ylabel('\DeltaF/F')
            xlabel('Time (seconds)')
            axis square;

        % Show tuning curve
        subplot(1,3,3); 
            plot(x,yMean,'-ok');
            hold on;
            errorbar(x(1:(end-1)),yMean(1:(end-1)),ySEM(1:(end-1)),'ok');
            for j = 1:(cells.uniqStims-1), plot(x(j),y(j,:),'o','Color',0.5*[1,1,1]); end
            set(gca,'Box','off','XTick',0:45:360);
            xlim([0,360]);
            title('Tuning curve');
            ylabel('\DeltaF/F')
            xlabel(sprintf('Stimulus direction (%s)', char(176)))
            axis square;    

        pause; drawnow;
        clf;
    end
end

%% Show maps of orientation and direction preference
if(showFunctionalMaps)
    meanCellResponses = mean(cells.meanStimResponse(:,1:(end-1),:),3);
    LUT           = hsv(360);
    markerMax     = 12;
    clippingLimit = 0.6; % clipping limit of selectivity (can range from 0 to 1).
    names         = {'Direction','Orientation'};
    figure('units','normalized','outerposition',[0 0 1 1])
    for isOrientation = 0:1
        % Compute orientation/direction preference using a vector sum approach
        theta       = wrapTo2Pi(angle(vectorSum(meanCellResponses,1+isOrientation)))*180/pi;
        offset      = (1+isOrientation);

        % Compute orientation/direction selectivity
        selectivity = abs(vectorSum(meanCellResponses,1+isOrientation))./sum(meanCellResponses,2);
        selectivity(selectivity>1 | selectivity<0) = 0;
        selectivity = selectivity/clippingLimit; selectivity(selectivity>1) = 1;

        % Show an orientation preference/direction map
        subplot(2,3,1+3*(1-isOrientation)); 
            imagesc(ones(512,512,3)); % First create a blank map to plot our cells on
            colormap(hsv); caxis([0 360/offset]); colorbar('Location','southoutside');
            title([names{1+isOrientation} ' preference map']); 
            axis image;
            hold on;
            for i = 1:cellNumber %% Now wel'll add each cell to our blank images
                markerSize = round(markerMax*selectivity(i));
                if(markerSize==0), continue; end
                plot(cells.xPos(i),cells.yPos(i),'ok','MarkerSize',markerSize,'MarkerFaceColor',LUT(1+floor(theta(i)) ,:));
            end

        % Show a histogram of tuning preference
        x = linspace(0,360/offset,9);
        subplot(2,3,2+3*(1-isOrientation)); 
            hBar = bar(x,hist(theta/offset,x)/length(theta));
            title('Histogram')
            ylabel('% cells');
            xlabel(sprintf('%s preference (\\circ)',names{1+isOrientation}));
            xlim([-22.5 (360+22.5)]/offset)
            axis square;
            set(gca,'Box','off');
            set(hBar,'BarWidth',0.8);

        % Show a histogram of tuning selectivity
        x = linspace(0,clippingLimit,9);
        subplot(2,3,3+3*(1-isOrientation)); 
            hBar =  bar(x,hist(selectivity,x)/length(selectivity));
            title('Histogram')
            ylabel('% cells');
            xlabel(sprintf('%s selectivity',names{1+isOrientation}));
            xlim([-0.05 1.05*clippingLimit])
            axis square;
            set(gca,'Box','off');
            set(hBar,'BarWidth',0.8);
    end
end