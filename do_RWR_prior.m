function [posteriors label_img] = do_RWR_prior(img,seeds,labels,c,lambda2,nei,sigma_c,isKeepConnect,imgName)
if nargin<6, nei = 0; end;
if nargin<7, sigma_c = 60; end;
if nargin<8, isKeepConnect = 1; end;
if nargin<9, imgName = 0; end;

[X Y Z]=size(img); N = X*Y; % image size
K = max(labels);            % number of labels

edgeRegion = ones(N,1);
%% find edge
% im = im2double(img);
% for i = 1:Z
%     grad(:,:,i)=abs(gradient(im(:,:,i)));
% end
% gradAvg = sum(grad,Z)./Z;
% edgeRegion = gradAvg>0.03;
% edgeRegion = imdilate(edgeRegion,ones(2,2));
% figure,imshow(edgeRegion);
% edgeRegion = reshape(edgeRegion,N,1);

%% Build graph with symmetric weight matrix W
if Z > 1,
    img = colorspace('Lab<-', img); % convert color space
end;
% img = double(img);
imgVals = reshape(img,N,Z);
[points edges] = lattice(X,Y,nei); clear points;
weights = makeweights(edges,imgVals,sigma_c);
edgeWeights = weights;
%% makeGMM
% addpath('temp');
if imgName == 0
%     fprintf('make GMM\n');
%     tic
    for k = 1:K
        seedsInd{k} =  seeds(labels==k);
    end
    w = makeGMMmutil(double(img),seedsInd);
%     toc
else if ~exist(['temp\',imgName,'_w.mat'])
        fprintf('make GMM\n');
        tic
        for k = 1:K
            seedsInd{k} =  seeds(labels==k);
        end
        w = makeGMMmutil(double(img),seedsInd);
        toc
        save(['temp\',imgName,'_w.mat'],'w');
    else
        load(['temp\',imgName,'_w.mat']);
    end
end
        
%% show GMM probablity
% Estimate likelihoods
likelihoods = zeros(N,K);
for k=1:K,
    likelihoods(:,k) = w(:,k);%/sum(w(:,k));
    w(:,k) = likelihoods(:,k);
end;

% % Estimate posteriors
% prob = sparse(1:N,1:N,1./sum(likelihoods,2))*likelihoods;
% [vals idx] = max(prob'); label_img = reshape(idx,X,Y);  clear vals idx;
% posteriors = reshape(prob,X,Y,K);
% figure; clf;set(gcf,'Position',[100,500,size(img,2)*(K+1),size(img,1)]);
% for k=1:K, 
%     prob_img = posteriors(:,:,k);
% %     prob_img = sc(posteriors(:,:,k),'prob_jet');
%     subplot(1,K+1,k); imagesc(prob_img); %clear prob_img;
% end;
% subplot(1,K+1,K+1); imshow(2-label_img);
% % imwrite(2-label_img,'temp.png');
%% initialization
W = adjacency(edges,weights,N); %clear edges weights;

    d = sum(W,2);
    D = diag(d);
    g= sum(w,2);
%     L = diag(d) - W;                         % Laplacian matrix
%     L = sparse(L);    
%     npt = N; 
    lines = zeros(N,K);
    
    for k=1:K,
        label_idx = find(labels(:)==k);
        Mk = size(label_idx,1);
        lines(seeds(label_idx(:)),k) = 1/Mk;
        clear label_idx;
        %test: choose connected region including seeds
        if isKeepConnect
%             edgeRegionT = reshape(edgeRegion,[X,Y]);
            tmpF = label_img==k;
            %tmpF =  tmpF | edgeRegionT;
%             figure; imshow(tmpF); 
            tmpL = bwlabel(tmpF,8);%label connected region
            LFR = tmpL(lines(:,k)>0);
            LFR = unique(LFR);
            if LFR(1)==0
                LFR(1)=[];
            end
            LFR_img = zeros(X,Y);
            for i = 1:length(LFR)
                LFR_img(tmpL==LFR(i))=1;
            end
            LFR_img = imdilate(LFR_img,ones(3));
%             figure; imshow(LFR_img);        
            edgeRegion = edgeRegion.*reshape(LFR_img,N,1);
        end
        
%         w(lines(:,k)>0,k)=g(lines(:,k)>0);
%         gk = w(:,k)+(g-w(:,k))/(K-1);
        gk = w; gk(:,k)=[]; gk = max(gk,[],2);
%         gk=g;
        Dg = sparse(1:N,1:N,gk.*edgeRegion);
        D_all = D+lambda2*Dg;
        A = D_all-(1-c)*W;
        b = (1-c)*lambda2*w(:,k).*edgeRegion+D_all*lines(:,k)*c;
%         tic
%         [VV,DD] = eigs(A,200,'sm');
%         R(:,k) = VV*(pinv(DD)*(VV'*b));  
%         toc
%         tic
        R(:,k) = A\b;
% %       R(:,k) = cs_lusol(A,b);
% 
%         toc
    end;
     
        %test
%     for k = 1:K
%         label_idx = find(labels(:)==k);
%         Mk = size(label_idx,1);
%         lines(seeds(label_idx(:)),k) = 1;
%         clear label_idx;
%     end
%     lines_all = sum(lines,2)>0;
%     for k = 1:K
%         Dg = sparse(1:N,1:N,g.*edgeRegion);
%         D_all = D+lambda2*Dg;
%         A = D_all-W;
%         b = lambda2*w(:,k).*(1-lines_all).*edgeRegion+lambda2*Dg*lines(:,k).*edgeRegion;
%         R(:,k) = A\b; 
%     end

%% output
% Estimate likelihoods
likelihoods = zeros(N,K);
for k=1:K,
    likelihoods(:,k) = R(:,k)/sum(R(:,k));
end;

% Estimate posteriors
prob = sparse(1:N,1:N,1./sum(likelihoods,2))*likelihoods;
[vals idx] = max(prob'); label_img = reshape(idx,X,Y);  clear vals idx;
posteriors = reshape(prob,X,Y,K);

end