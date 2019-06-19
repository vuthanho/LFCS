close all
folder = strcat(pwd,'/../LFColorSample/EXR/images/');

im{1} = double(exrread(strcat(folder,'028.exr')));
im{2} = double(exrread(strcat(folder,'027.exr')));
im{3} = double(exrread(strcat(folder,'019.exr')));
im{4} = double(exrread(strcat(folder,'011.exr')));
im{5} = double(exrread(strcat(folder,'018.exr')));
im{6} = double(exrread(strcat(folder,'010.exr')));

values = cell(1,6);
factor = [1 .25];
output_size = size(imresize(im{1},factor(1),'nearest') );

for i=1:6
    values{i}.I1exp0   = zeros(output_size);
    values{i}.H        = zeros(3, 3);
    values{i}.gamma    = zeros(1, 2);
    values{i}.max_clip = zeros(1, 3);
    values{i}.min_clip = zeros(1, 3);
end

gamma = 1;
coef = [];
clip = [0 0.99];
use_sift = [1 0];

values{2} = compute_colorstabilization( double( imresize(im{2}, factor(1),'nearest', 'Colormap', 'original') ), ...
                double( imresize(im{2}, factor(2),'nearest', 'Colormap', 'original') ),  ...
                double( imresize(im{1}, factor(2),'nearest', 'Colormap', 'original') ),  ...
                coef, clip, gamma, use_sift, [] );

%% first column

values{5} = compute_colorstabilization( double( imresize(im{5}, factor(1),'nearest', 'Colormap', 'original') ), ...
                double( imresize(im{5}, factor(2),'nearest', 'Colormap', 'original') ),  ...
                double( imresize(values{2}.I1exp0,factor(2) )),  ...
                coef, clip, gamma, use_sift, [] );

values{6} = compute_colorstabilization( double( imresize(im{6}, factor(1),'nearest', 'Colormap', 'original') ), ...
                double( imresize(im{6}, factor(2),'nearest', 'Colormap', 'original') ),  ...
                double( imresize(values{5}.I1exp0, factor(2),'nearest', 'Colormap', 'original') ),  ...
                coef, clip, gamma, use_sift, [] );            

%% second column

values{3} = compute_colorstabilization( double( imresize(im{3}, factor(1),'nearest', 'Colormap', 'original') ), ...
                double( imresize(im{3}, factor(2),'nearest', 'Colormap', 'original') ),  ...
                double( imresize(values{2}.I1exp0, factor(2),'nearest', 'Colormap', 'original') ),  ...
                coef, clip, gamma, use_sift, [] );

values{4} = compute_colorstabilization( double( imresize(im{4}, factor(1),'nearest', 'Colormap', 'original') ), ...
                double( imresize(im{4}, factor(2),'nearest', 'Colormap', 'original') ),  ...
                double( imresize(values{3}.I1exp0, factor(2),'nearest', 'Colormap', 'original') ),  ...
                coef, clip, gamma, use_sift, [] );

im4nobright = reshape((values{4}.H*(reshape(double( imresize(im{4},factor(1),'nearest' , 'Colormap', 'original') ),[],3)...
                .^(values{4}.gamma(1) / values{4}.gamma(2)))')',size(values{4}.I1exp0));

cmp01 = values{3}.I1exp0;
cmp01(:,1:1920/2,:) = values{5}.I1exp0(:,1:1920/2,:);            
            
cmp1 = values{4}.I1exp0;
cmp1(:,1:1920/2,:) = values{6}.I1exp0(:,1:1920/2,:);

cmp2 = im4nobright;
cmp2(:,1:1920/2,:) = values{6}.I1exp0(:,1:1920/2,:);

imshow(cmp01)
title('5-3')
figure
imshow(cmp1)
title('6-4')
% figure
% imshow(cmp2)
% title('6-4 (no intensity consideration)')