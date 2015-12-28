clear ;
close all;

addpath 'algorithms'
out = ['results\'];
if ~exist(out)
    mkdir(out);
end

%% parameters            
only_name= '41004';
img_name = ['./imgs/' only_name '.jpg'];
ref_name = ['./scribbles/' only_name '.bmp'];


nei = 1;            % 0: 4-neighbors, 1: 8-neighbors
c = 0.0004;%1e-3;%         % restarting probability of RWR
sigma_c = 60;       % color variance
scale = 1.0;        % image resize
lambda = 2e-10;     % parameter for unitary
isKeepConnect = 0;  % 1: only consinder the connected regions with seeds;  0: otherwise.
reset(RandStream.getGlobalStream); % fix the random seed for initalization of GMM

saveProb = 0; % 1: save the probability image
%% main routine
img = imread(img_name); img = imresize(img,scale);
[K, labels, idx] = seed_generation(ref_name,scale);

%% RWR with prior
% run('vlfeat-0.9.13/toolbox/vl_setup');
st=clock;
[posteriors label_img] = do_RWR_prior(img,idx,labels,c,lambda,nei,sigma_c,isKeepConnect);
fprintf('subRW took %.2f second\n',etime(clock,st));

% display
[imgMasks,segOutline,imgMarkup]=segoutput(im2double(img),label_img); %clear imgMasks segOutline;

outPath = [out,'ours\'];
if ~exist(outPath)
    mkdir(outPath);
end

figure; clf;set(gcf,'Position',[100,500,size(img,2)*(K+1),size(img,1)]);
for k=1:K 
    prob_img = sc(posteriors(:,:,k),'prob_jet');
    if saveProb == 1
        imwrite(prob_img,[outPath,only_name,'_prob',num2str(k),'.png']);
    end
    subplot(1,K+1,k); imshow(prob_img); clear prob_img;
end;
subplot(1,K+1,K+1); imshow(imgMarkup);
figure,imshow((K-imgMasks)/(K-1));

imwrite(imgMarkup,[outPath,only_name,'_bound.png']);
imwrite((K-imgMasks)/(K-1),[outPath,only_name,'_binary.png']);