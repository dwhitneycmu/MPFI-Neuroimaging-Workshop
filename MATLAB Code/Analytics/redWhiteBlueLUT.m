function LUT= redWhiteBlueLUT(n)
% Create red-white-blue colormap

if(nargin<1), n = 64; end

% Define start, middle, end colors
r = [1 0 0]; %# start
w = [1 1 1]; %# middle
b = [0 0 1]; %# end

%# colormap of size 64-by-3, ranging from red -> white -> blue
nUnits = round(n/2);
c1 = zeros(nUnits,3); c2 = zeros(nUnits,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), nUnits);
    c2(:,i) = linspace(w(i), b(i), nUnits);
end
LUT = flipud([c1(1:end-1,:);c2]);