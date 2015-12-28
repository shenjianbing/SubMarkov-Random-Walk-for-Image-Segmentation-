function [posteriors label_img] = do_segmentation(img,seeds,labels,c,nei,sigma_c)

if nargin<5, nei = 0; end;
if nargin<6, sigma_c = 60; end;

[X Y Z]=size(img); N = X*Y; % image size
K = max(labels);            % number of labels

% Build graph with symmetric weight matrix W
if Z > 1,
    img = colorspace('Lab<-', img); % convert color space
end;
imgVals = reshape(img,N,Z);
[points edges] = lattice(X,Y,nei); clear points;
weights = makeweights(edges,imgVals,sigma_c);
W = adjacency(edges,weights,N); clear edges weights;

% Random Walks with Restart (RWR)
E = sparse(1:N,1:N,ones(N,1)); iD = sparse(1:N,1:N,1./sum(W));
lines = zeros(N,K);
for k=1:K,
    label_idx = find(labels(:)==k);
    Mk = size(label_idx,1);
    lines(seeds(label_idx(:)),k) = 1\Mk;
    clear label_idx;
end;

P = iD*W; 
R = c*(E-(1-c)*P)\lines;

% Estimate likelihoods
likelihoods = zeros(N,K);
for k=1:K,
    likelihoods(:,k) = R(:,k)/sum(R(:,k));
end;

% Estimate posteriors
prob = sparse(1:N,1:N,1./sum(likelihoods,2))*likelihoods;
[vals idx] = max(prob'); label_img = reshape(idx,X,Y);  clear vals idx;
posteriors = reshape(prob,X,Y,K);

