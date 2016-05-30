function [permutation, Idata, rbfWeights, rbfCenters] = ...
    pmComputeSubSet(Adata, ...            % start photon position 2xN
                    Bdata, ...            % end photon position 2xN
                    subset_size, ...      % size of subset in subset assignment
                    beta)                 % weight between aubset assignment energy terms
%PMCOMPUTESUBSET Computes a matching between two point sets. Also returns
%the RBF weights and the centers

    %% error checking
    if size(Adata) ~= size(Bdata)
        error('Sizes of Adata and Bdata do not match in pmComputeMatching.');
    end

    %% Default values
    if nargin < 3
        subset_size=300;
    end
    if nargin < 4
        beta=0.005;
    end
    
    %% pick subset for assignment problem
    N = size(Adata,2);     % number of samples
    M = min(subset_size, N);
    Asub = Adata(:,1:M);   % samples are in random order, so pick first M
    Bsub = Bdata(:,1:M);
    c = pmComputeFullSet(Asub, Bsub, beta); % compute assignment
    Bsub = Bsub(:,c);      % reorder Bsub to have the indices match
    
    %% setup the rbf system
    V = pmRBF(Asub', Asub');
    rbfWeights=V\Bsub';    % compute weights
    rbfCenters = Asub;
    
    %% compute target positions from rbf
    Vall = pmRBF(rbfCenters', Adata');
    Idata = Vall * rbfWeights; % morph full point set
    Idata = Idata';
    
    %% do a greedy refinement of the assignment
    permutation = pmRefinement(Idata, Bdata);

end

