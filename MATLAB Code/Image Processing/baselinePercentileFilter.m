function highpassFilteredTrace = baselinePercentileFilter(inputTrace,fps,filteredCutOff,desiredPercentileRank)
% highpassFilteredTrace = baselinePercentileFilter(inputTrace,fps,filteredCutOff,desiredPercentileRank)
%
% Uses a combination of 1-d median/butterworth high-pass filtered to compute a baseline 
% for the input trace:
% *fps - sampling rate (frames per second)
% *filteredCutOff - low-pass filter cutoff (in seconds)
% *filterType - which filtering algorithm to use (fastMedFilt1d or rankOrder). 
%               Default filter is fastMedFilt1d (which is much faster).
% desiredPercentileRank - if filter is rankOrder, then the percentile cutoff
%                         can be picked (default is 50% for computing median).
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2016.

if(nargin<2), fps = 15.1515;                end % in Hz
if(nargin<3), filteredCutOff = 60;          end % in seconds
if(nargin<4), desiredPercentileRank = 50;   end % only valid if rankOrder
isVerbose = false;

% Compute a low-pass median filter
paddingLength = ceil(length(inputTrace)/1);
paddedTrace   = [inputTrace(paddingLength:-1:1); inputTrace; inputTrace(paddingLength:-1:1)];
%filteredTrace = rankOrderFilter(paddedTrace,round(filteredCutOff*fps),desiredPercentileRank);
filteredTrace = percentileFilt1(paddedTrace,desiredPercentileRank,round(filteredCutOff*fps));
filteredTrace = filteredTrace(paddingLength+[1:length(inputTrace)]);
if(isVerbose), figure; plot(filteredTrace); end

% The low-pass filter the filtered trace to smooth it out
butterWorthOrder = 1;
Wn = (1/filteredCutOff)/(fps/2);% normalized cutoff frequency in unit of nyquist frequency
[b,a] = butter(butterWorthOrder, Wn, 'low');
highpassFilteredTrace = filtfilt(b,a,[filteredTrace(paddingLength:-1:1); filteredTrace; filteredTrace(1:paddingLength)]);
highpassFilteredTrace = highpassFilteredTrace(paddingLength+[1:length(inputTrace)]);
if(isVerbose), figure; plot(highpassFilteredTrace); end

function y = rankOrderFilter(x, N, p)
%RankOrderFilter Rank order filter for 1D signals
%  y = RankOrderFilter(x, window, thd) runs a rank-order filtering of order
%  N on x. y is the same size as x. To avoid edge effects, the x is expanded
%  by repeating the first and the last samples N/2 times. if x is a matrix,
%  RankOrderFilter operates along the columns of x.
%
%  Rank-order filter calculates the p'th percentile of the data on a N
%  sized window round each point of x. p can be a number between 0 and 100.
%
%  When p is equal to 50, the output of this function will be the same as 
%  MATLAB's PRCTILEFILT1(x,N); however, RankOrderFilter is almost always much
%  faster and needs less memory. 
%
%  When p is close to 0 (or to 100), a RankOrderFilter calculates an
%  approximate lower (or upper) envlope of the signal.
%
%  Copyright 2008, Arash Salarian
%  mailto://arash.salarian@ieee.org
%

if isrow(x);
	x=x';
	transposedx=true;
else 
	transposedx=false;
end

[m, n] = size(x);
y = zeros(m,n);
if rem(N,2) == 1
    k = (N-1)/2;
else
    k = N /2;
end

for i=1:n
    X = [x(1,i)*ones(k,1);x(:,i);x(end,i)*ones(k,1)];
     for j=1:m
         y(j,i) = percentile(X(j:j+N-1), p);
    end
end
if transposedx
	y=y';
end

% Percentile calculated the k'th percentile of x. This function is similar 
% to, but generally much faster than MATLAB's prctile function.
function y = percentile(x, k)
x = sort(x);
n = size(x,1);

p = 1 + (n-1) * k / 100;

if p == fix(p)
    y = x(p);
else
    r1 = floor(p); r2 = r1+1;
    y = x(r1) + (x(r2)-x(r1)) * k / 100;
end

function y = percentileFilt1(x,percentile,n,blksz,DIM)
%percentileFilt1  One dimensional median filter.
%   Y = percentileFilt1(X,percentile,N) returns the output of the order N, one dimensional
%   percentile filtering of X.  Y is the same size as X; for the edge points,
%   zeros are assumed to the left and right of X.  If X is a matrix,
%   then percentileFilt1 operates along the columns of X.
%
%   If you do not specify N, percentileFilt1 uses a default of N = 3.
%   For N odd, Y(k) is the percentile of X( k-(N-1)/2 : k+(N-1)/2 ).
%   For N even, Y(k) is the percentile of X( k-N/2 : k+N/2-1 ).
%
%   Y = percentileFilt1(X,N,BLKSZ) uses a for-loop to compute BLKSZ ("block size") 
%   output samples at a time.  Use this option with BLKSZ << LENGTH(X) if 
%   you are low on memory (percentileFilt1 uses a working matrix of size
%   N x BLKSZ).  By default, BLKSZ == LENGTH(X); this is the fastest
%   execution if you have the memory for it.
%
%   For matrices and N-D arrays, Y = percentileFilt1(X,N,[],DIM) or 
%   Y = percentileFilt1(X,N,BLKSZ,DIM) operates along the dimension DIM.
%   
%   % Example:
%   %   Construct a noisy signal and apply a 10th order one-dimensional 
%   %   median filter to it.
%
%   fs = 100;                               % Sampling rate                                   
%   t = 0:1/fs:1;                           % Time vector
%   x = sin(2*pi*t*3)+.25*sin(2*pi*t*40);   % Noise Signal - Input
%   y = percentileFilt1(x,50,10);           % Median filtering - Output
%   plot(t,x,'k',t,y,'r'); grid;            % Plot 
%   legend('Original Signal','Filtered Signal')
%
%   See also MEDIAN, FILTER, SGOLAYFILT, and MEDFILT1 in the Image
%   Processing Toolbox.

%   Author(s): L. Shure and T. Krauss, 8-3-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.6 $  $Date: 2012/10/29 19:31:41 $

% Validate number of input arguments
narginchk(1,5);
if nargin < 2, percentile = 50; end
if nargin < 3, n = []; end
if nargin < 4, blksz = []; end
if nargin < 5, DIM = []; end

% Check the input data type. Single precision is not supported.
% try
%     chkinputdatatype(x,n,blksz,DIM);
% catch ME
%     throwAsCaller(ME);
% end

% Check if the input arguments are valid
if isempty(n)
  n = 3;
end

if ~isempty(DIM) && DIM > ndims(x)
	error(message('signal:medfilt1:InvalidDimensions'))
end

% Reshape x into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first (row) dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

% Verify that the block size is valid.
siz = size(x);
if isempty(blksz),
	blksz = siz(1); % siz(1) is the number of rows of x (default)
else
	blksz = blksz(:);
end

% Initialize y with the correct dimension
y = zeros(siz); 

% Call medfilt1D (vector)
for i = 1:prod(siz(2:end)),
	y(:,i) = prctilefilt1d(x(:,i),n,blksz,percentile);
end

% Convert y to the original shape of x
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end


%-------------------------------------------------------------------
%                       Local Function
%-------------------------------------------------------------------
function y = prctilefilt1d(x,n,blksz,percentile)
%PRCTILEFILT1D  One dimensional median filter.
%
% Inputs:
%   x     - vector
%   n     - order of the filter
%   blksz - block size

nx = length(x);
if rem(n,2)~=1    % n even
    m = n/2;
else
    m = (n-1)/2;
end
X = [zeros(m,1); x; zeros(m,1)];
y = zeros(nx,1);

% Work in chunks to save memory
indr = (0:n-1)';
indc = 1:nx;
for i=1:blksz:nx
    ind = indc(ones(1,n),i:min(i+blksz-1,nx)) + ...
          indr(:,ones(1,min(i+blksz-1,nx)-i+1));
    xx = reshape(X(ind),n,min(i+blksz-1,nx)-i+1);
    y(i:min(i+blksz-1,nx)) = prctile(xx,percentile,1);
end
