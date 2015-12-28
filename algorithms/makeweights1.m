function [edgeWeights,nullWeights,edgeFeatures, edgeNodeIndex]=makeweights1(img)


A = im2double(img);
N = size(A, 1);
M = size(A, 2);

% make set of pixels of the image
Av = reshape(A, [], 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create edge classes

% We start by taking differences of adjacent pixels. The order here is the
% same that is used in the C++ code provided online, it is important that
% the order of the entries is retained.
% This is done in four blocks (directions).

nodeId = reshape(1 : N * M, [N, M]);

% vertical
% create elist-part and use that to compute gradients
tmp1 = nodeId(1 : end - 1, :);
tmp2 = nodeId(2 : end, :);
edgeNodeIndex = [ tmp1(:), tmp2(:) ];

% horizontal, to right and to left
tmp1 = nodeId(:, 2 : end);
tmp2 = nodeId(:, 1 : end - 1);
edgeNodeIndex = [ edgeNodeIndex; tmp1(:), tmp2(:) ];

% diagonal \
tmp1 = nodeId(1 : end - 1, 1 : end - 1);
tmp2 = nodeId(2 : end, 2 : end);
edgeNodeIndex = [ edgeNodeIndex; tmp1(:), tmp2(:) ];

% diagonal /
tmp1 = nodeId(2 : end, 1 : end - 1);
tmp2 = nodeId(1 : end - 1, 2 : end);
edgeNodeIndex = [ edgeNodeIndex; tmp1(:), tmp2(:) ];

% make feature vectors for each edges
edgeFeatures = Av(edgeNodeIndex(:, 1), :) - Av(edgeNodeIndex(:, 2), :);

% clear edgeNodeIndex

% now we are done creating the vectors of pixel differences
% next we transform it to the standard (Gaussian) edge weights
edgeWeights = sum( edgeFeatures .^ 2, 2);
nullWeights = (edgeWeights < eps);

sigma = mean( edgeWeights );
edgeWeights = 0.05 + 0.95 * exp( -edgeWeights / (2 * sigma));




