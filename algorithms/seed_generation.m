function [nlabels, labels, idx] = seed_generation(ref_name,scale)

if nargin<2, scale = 1; end;

ref=im2double(imread(ref_name)); ref = imresize(ref,scale);
L{1} = find(ref(:,:,1)==1.0 & ref(:,:,2)==0.0 & ref(:,:,3)==0.0); % R
L{2} = find(ref(:,:,1)==0.0 & ref(:,:,2)==1.0 & ref(:,:,3)==0.0); % G
L{3} = find(ref(:,:,1)==0.0 & ref(:,:,2)==0.0 & ref(:,:,3)==1.0); % B

num = 0;
nlabels = 0;
for i=1:3
    nL = size(L{i},1);
    if nL > 0
        nlabels = nlabels + 1;
        labels(num+1:num+nL) = nlabels;
        idx(num+1:num+nL) = L{i};
        num = num + nL;
    end;
end;
