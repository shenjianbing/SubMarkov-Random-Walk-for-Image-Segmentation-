function [w] = makeGMMmutil(A, seeds)

if ~exist('EM.m','file')
    addpath('GMM-GMR-v2.0');
end
% addpath('../em/netlib');

A = reshape(A, [], 3)';
A = A;%+0.00001*randn(size(A));

ncenters = 5;
K = length(seeds);

n = size(A,2);
stabF = 0.01;
for k = 1:K
    fgPix = seeds{k};   

    X = A(:,fgPix);
    [Unaries, Mu, Sigma] = EM_init_kmeans(X, ncenters);
    [unarypotF, muF, sigmaF] = EM(X, Unaries, Mu, Sigma);

    

    wb = zeros(n,ncenters);   

    for i=1:ncenters        

        L = chol(inv(sigmaF(:,:,i)+stabF*eye(3)));
        tmp = L*(A - repmat(muF(:,i),1,n));
        wb(:,i) = -log(unarypotF(i)) + 0.5*log(det(sigmaF(:,:,i) + stabF*eye(3))) + 0.5* sum( tmp.^2, 1)'; 
    end
    
    wb2 = min(wb, [], 2) + 3/2 * log(2*pi);

    
%     while min(wb2) < 0
%         stabF = stabF + 0.1;
%         for i=1:ncenters
%             L = chol(inv(sigmaF(:,:,i)+stabF*eye(3)));
%             tmp = L*(A - repmat(muF(:,i),1,n));
%             wb(:,i) = -log(unarypotF(i)) + 0.5*log(det(sigmaF(:,:,i) + stabF*eye(3))) + 0.5* sum( tmp.^2, 1)';
%         end
%         wb2 = min(wb, [], 2) + 3/2 * log(2*pi);        
%     end    
    w(:,k) = max(100-wb2,eps);
%     w(:,k) = 100-wb2;
%     w(fgPix,k) = 100;
end