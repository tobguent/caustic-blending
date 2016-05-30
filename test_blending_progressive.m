%% Description
% This file is an extension of "test_blending.m". Here, we also sample
% two images and compute a blending, but here, this is done progressively
% by iteratively generating further point sets. The density estimates use
% progressively shrinking filter kernels, which results in unbiased
% solutions in the first and last frame, i.e., the start and end frame a
% perfectly interpolated. This script exports a file series.
% In the default setting the number of progressive iterations is
% 'progressive_its = 100'. Depending on your machine, the calculation
% might take minutes (5 minutes on our test machine for 100 iterations).
% Increase the number of progressive iterations to obtain less noisy
% results.

%% clean up
clear;
clf;

%% compile the mex files
pmCompileMexFiles;

%% experiment parameters
fullset_size = 20000;   % total size of the point sets
subset_size = 400;      % subset size for partial matching
beta = 0.04;            % energy weight between 'short path' and 'structure preservation'
ratweight = 0.2;        % blending weight that makes the paths more turbulent
progressive_its = 100; % number of progressive iterations

radiusSquared = 6*6;    % initial radius
M = 500;  N = 500;      % resolution of the resulting animation
num_frames = 120;       % number of frames in the animation
densities = zeros(M,N,num_frames);  % output buffer that stores the density estimates for each frame of the animation

%% create CDF for source distribution
imageA = imread('images/caustic.png');
imageA = mean(imageA,3) / 255;
[colCDFA, rowCDFA, probA] = pmCreateCDF2d(imageA);   % can be done once

%% create CDF for target distribution
imageB = imread('images/flower.png');
imageB = mean(imageB,3) / 255;
[colCDFB, rowCDFB, probB] = pmCreateCDF2d(imageB);   % can be done once

tic;
for pro = 1 : progressive_its
   
    % generate a point set for the source
    rnd = rand(4, fullset_size);     % generate 4 random numbers in [0,1] for each sample
    [Adata, Aflux] = pmSampleInvCDF2d(colCDFA, rowCDFA, probA, imageA, rnd);
    
    % generate a point set for the target
    rnd = rand(4, fullset_size);     % generate 4 random numbers in [0,1] for each sample
    [Bdata, Bflux] = pmSampleInvCDF2d(colCDFB, rowCDFB, probB, imageB, rnd);

    if pro == 1
       % first iteration constructs the RBF weights
       [permutation, Idata, rbfWeights, rbfCenters] = pmComputeSubSet(Adata, Bdata, subset_size, beta); 
    else
        % later iterations use the precomputed RBF weights
        Vall = pmRBF(rbfCenters', Adata');
        Idata = Vall * rbfWeights; % interpolate full point set
        Idata = Idata';
        % do a greedy refinement of the assignment
        permutation = pmRefinement(Idata, Bdata);
    end
    
    % calculate the density frame for each frame and sum results up
    for it = 1:num_frames
        blend = it / num_frames;
        [Tdata, Tflux] = pmMorph(Adata, Idata, Bdata, permutation, Aflux, Bflux, blend, ratweight);
        densities(:,:,it) = densities(:,:,it) + pmDensityEstimate( Tdata, Tflux, sqrt(radiusSquared), M, N) ./ progressive_its;
        
        % shrink the radius of the density estimator (according to Knaus-Zwicker)
        radiusSquared = radiusSquared * (it + 0.6666666) / (it + 1);
    end
    
end
toc;

%% Display source and target image
subplot(1,3,1);
imshow(imageA);
title('Source image');

subplot(1,3,2);
imshow(imageB);
title('Target image');


%% Make a video of the blended image
for it = 1:num_frames
    blend = it / num_frames;
    
    subplot(1,3,3);
    imshow(densities(:,:,it));
    title('Blended density');
    
    % export the image
    imwrite(densities(:,:,it), sprintf('sequence%03d.bmp',it));
    
    pause(0.0033); % small delay between consecutive frames
end
