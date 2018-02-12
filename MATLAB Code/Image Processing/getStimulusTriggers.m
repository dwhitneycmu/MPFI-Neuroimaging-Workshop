function [stimOnsetTime,stimID,frameTimes]=getStimulusTriggers(path)
% [stimOnsetTime,stimID,meanCCDTime]=getStimulusTriggers(path)
%
% Helper function to load timing information about the image stack relative
% to trials blocks of the visual stimulus. Functions take in the path of
% the triggers (Spike2), and returns:
% 1.) stimOnsetTime - When a visual stimulus was presented (in seconds)
% 2.) stimID        - The identity of the stimulus
% 3.) frameTimes    - Time (in seconds) when each imaging frame was
%                     acquired (relative to stimOnsetTime)
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

% Read file with frame times (if present)
triggerPath = [path 'frametrigger.txt'];
if(exist(triggerPath) ~= 0)  
    CCD_FramesTimes  = load(triggerPath); 
    startTime        = CCD_FramesTimes(1);
    frameTimes       = CCD_FramesTimes-startTime;
else
    startTime  = [];
    frameTimes = [];
end

% Read file with frame times for stimulus onset (and stimulus label)
stimulusOnPath = [path 'stimontimes.txt'];
if(exist(stimulusOnPath) ~= 0)
    stimData = load(stimulusOnPath);
    if(~isempty(stimData))
        if stimData(1)==0
            stimID        = stimData(3:2:end);
            stimOnsetTime = stimData(4:2:end);
        else
            stimID        = stimData(1:2:end);
            stimOnsetTime = stimData(2:2:end);
        end
        stimOnsetTime = stimOnsetTime-startTime;
    end
else
    stimID      = [];
    stimOnsetTime  = [];
end
return
