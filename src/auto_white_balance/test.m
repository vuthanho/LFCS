%% main to test awb.m

clear all
close all

folder      = 'D:\Cours\2A\Internship\LFcolorStable_olivier\LF_color_sample_data\JPG\images';
im_shifted  = im2double(imread(strcat(folder,'\027.jpg')));
im          = im2double(imread(strcat(folder,'\028.jpg')));
factor      = [1 1];
gamma       = 2.3142;
use_sift    = [1 0];
clip        = [0 1];


im_med          = imresize(im, factor(2), 'Colormap', 'original');
im_shifted_med  = imresize(im_shifted, factor(2), 'Colormap', 'original');
I1tmp = im_shifted_med;
I2tmp = im_med;
% Correspondences points computing
[pos1, pos2, I1_reshape, I2_reshape] = find_correspondencesSIFT( (I1tmp./max(I1tmp(:))), (I2tmp./max(I2tmp(:))), I1tmp, I2tmp, use_sift(2));

subplot (1,2,1);
imshow (I2tmp);
hold on;
plot (pos2(2,:), pos2(1,:), 'r*');
title('Ref')

subplot (1,2,2);
imshow (I1tmp);
hold on;
plot (pos1(2,:), pos1(1,:), 'b*');
title('Test')
drawnow

value_27 = compute_colorstabilization( im_shifted_med, ...
                im_shifted_med, im_med,  ...
                [], clip, gamma, use_sift );

value_28 = compute_colorstabilization( im_med, ...
                im_med, im_med,  ...
                [], clip, gamma, use_sift );
            
im_shifted_CS   = value_27.I1exp0.^(1/value_28.gamma(2));
im_CS           = value_28.I1exp0.^(1/value_28.gamma(2));

figure
subplot(121)
imshow(im_CS)
title('Ref')
subplot(122)
imshow(im_shifted_CS)
title('Test')