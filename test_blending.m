%% Description
% This example samples two images with an equal number of points and 
% displays a blending between the two sets via animation. The resulting 
% point sets are rasterized by a density estimate for display.
% Note that only a finite number of samples is used. Thus, the result is
% very noisy. See "test_blending_progressive.m" for a progressive version
% that computes in the limit a noise free solution.

%% clean up
clear;
clf;

%% compile the mex files
pmCompileMexFiles;

%% experiment parameters
fullset_size = 50000;   % total size of the point sets
subset_size = 400;      % subset size for partial matching
beta = 0.04;            % energy weight between 'short path' and 'structure preservation'
ratweight = 0.2;        % blending weight that makes the paths more turbulent

%% create CDF for source distribution and sample point set for the source
imageA = imread('images/caustic.png');
imageA = mean(imageA,3) / 255;
[colCDFA, rowCDFA, probA] = pmCreateCDF2d(imageA);   % can be done once
rnd = rand(4, fullset_size);     % generate 4 random numbers in [0,1] for each sample
[Adata, Aflux] = pmSampleInvCDF2d(colCDFA, rowCDFA, probA, imageA, rnd);

%% create CDF for target distribution and sample point set for the target
imageB = imread('images/flower.png');
imageB = mean(imageB,3) / 255;
[colCDFB, rowCDFB, probB] = pmCreateCDF2d(imageB);   % can be done once
rnd = rand(4, fullset_size);     % generate 4 random numbers in [0,1] for each sample
[Bdata, Bflux] = pmSampleInvCDF2d(colCDFB, rowCDFB, probB, imageB, rnd);

%% Compute the matching and determine the RBF coefficients
[permutation, Idata, rbfWeights, rbfCenters] = pmComputeSubSet(Adata, Bdata, subset_size, beta);

%% output the point sets of source and target
subplot(2,3,1);
plot(Adata(1,:), Adata(2,:), '.r');
axis([-1 1 -1 1])
axis equal;
title('Source points: A');

subplot(2,3,2);
plot(Bdata(1,:), Bdata(2,:), '.r');
axis([-1 1 -1 1])
axis equal;
title('Target points: B');

%% output a density estimate of before and after
radius = 2.0;
M = 500;  N = 500;

subplot(2,3,4);
imshow(pmDensityEstimate( Adata, Aflux, radius, M, N));
title('Source density');

subplot(2,3,5);
imshow(pmDensityEstimate( Bdata, Bflux, radius, M, N));
title('Target density');

%% Make a video
num_frames = 120;
for it = 1:num_frames
    
    %% blend the photon positions and flux for the given time step
    blend = it / num_frames;
    [Tdata, Tflux] = pmMorph(Adata, Idata, Bdata, permutation, Aflux, Bflux, blend, ratweight);
    
    %% plot both the points, as well as a density estimate for the frame
    subplot(2,3,3);
    plot(Tdata(1,:), Tdata(2,:), '.r');
    axis([-1 1 -1 1])
    axis equal;
    title('Blended points');
    
    subplot(2,3,6);
    imshow(pmDensityEstimate( Tdata, Tflux, radius, M, N));
    title('Blended density');
    
    pause(0.0033);  % small delay between consecutive frames
end
