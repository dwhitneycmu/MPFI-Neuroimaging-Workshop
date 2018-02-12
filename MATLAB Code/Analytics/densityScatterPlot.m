function [h,histogramValues] = densityScatterPlot(x,y,xBins,yBins,isVerbose)
% function [h,histogramValues] = densityScatterPlot(x,y,xBins,yBins)
% Compute Density Scatter Plot on an Image
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2013.

if(nargin<5)
    isVerbose = true;
end
if(nargin<3)
    xBins = linspace(min(x(:)),max(y(:)),1000);
    yBins = linspace(min(y(:)),max(y(:)),1000);
end
edges{1} = xBins;
edges{2} = yBins;
histogramValues = hist3([x,y],'Edges',edges);

if(isVerbose)
    h = surf(xBins,yBins,histogramValues','EdgeColor','none'); view(2); axis square; colorbar;
    xlim([min(xBins) max(xBins)]); ylim([min(yBins) max(yBins)]);
else
    h = [];
end