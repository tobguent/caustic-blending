function [ D ] = pmDensityEstimate( pnts, flux, kernel, M, N )
%PMDENSITYESTIMATE Density estimation on MxN grid. Assumes domain size: [-1,1]x[-1,1]

% if not specified otherwise, use quadratic density map
if nargin < 5
    N = M;
end

% get coordinates
pnts = [0,-1; 1,0] * pnts; % rotate to have the axis matching on output
x = pnts(1,:);
y = pnts(2,:);

minx = -1;
miny = -1;
maxx = 1;
maxy = 1;

% compute splatting index
rx = int32((x - minx)/(maxx-minx) * (M))+1;
ry = int32((y - miny)/(maxy-miny) * (N))+1;

% size of an area element:
A = 2/M * 2/N;

% splat flux into the density map via scattering operation
D = zeros(M,N);
D = pmScatter(D, rx, ry, flux) ./ A;

% post-convolve with a kernel
if kernel > 0
    G = fspecial('disk', kernel);
    D = imfilter(D,G,'same');
end


end

