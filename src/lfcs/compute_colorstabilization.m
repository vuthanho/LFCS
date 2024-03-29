function values = compute_colorstabilization( I1, I1tmp, I2tmp, clip_v, gamma_ref, use_sift, ref_ref, limit )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Color matrix and gamma estimation given two image filename
%% I1 --> I2 (I2 is set as the reference)
%%
%% Inputs:  1. I1 -> char/matrix source original full size image
%%          2. I1tmp -> char/matrix source 1/4 smoothed image
%%          3. I2tmp -> char/matrix reference 1/4 smoothed image
%%          4. coef -> struct for ransac
%%          5. clip_v -> 1x2 array for min and max values
%%          6. gamma_ref -> gamma reference value
%%          7. use_sift -> 1-compute sift correspondences / 0-whole image pixels
%%          8. ref_ref -> gamma and H lastly used by ref
%%          9. limit -> limit after which the homography will be less applied
%%
%% Outputs: 1. values -> struct with -H (homography), 
%%                                    -gammas (estimated gamma values) 
%%                                    -I12 (linearized I1 image)
%%                                    -mask (binary image 1-pixels used for 
%%                                           color stabilization)
%%                                    -clipping (new clipping range min and max)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If SIFT, calculate correspondences using SIFT, otherwise consider 
%% all pixels in the image as one-to-one correspondences
if use_sift(1)
    [pos1, pos2, I1_reshape, I2_reshape] = find_correspondencesSIFT( (I1tmp./max(I1tmp(:))), (I2tmp./max(I2tmp(:))), I1tmp, I2tmp, use_sift(2));
else
    I1_reshape = reshape(I1tmp, [], 3);
    I2_reshape = reshape(I2tmp, [], 3);
end
% D = sqrt( (pos1(1,:)-pos2(1,:)).^2 + (pos1(2,:)-pos2(2,:)).^2 );
% index = D>std(D);
% I1_reshape = I1_reshape(~index,:);
% I2_reshape = I2_reshape(~index,:);
% pos1 = pos1(:,~index);
% pos2 = pos2(:,~index);



%% Clip values outside the given range
[r_l1, r_u1] = discard_saturated( I1_reshape, clip_v );
[r_l2, r_u2] = discard_saturated( I2_reshape, clip_v );

I1_reshape( r_u1 | r_l1 | r_u2 | r_l2, : ) = [ ]; 
I2_reshape( r_u1 | r_l1 | r_u2 | r_l2, : ) = [ ]; 
if use_sift(1)
    pos1(:, r_u1 | r_l1 | r_u2 | r_l2 ) = [ ]; 
    pos2(:, r_u1 | r_l1 | r_u2 | r_l2 ) = [ ]; 
end

[I1_reshape, ~] = clustering( I1_reshape, 2., -1e8 );
[I2_reshape, ~] = clustering( I2_reshape, 2., -1e8 );

count = 1;
for i=1:size(I1_reshape, 1)
    if( I1_reshape(i, 1) == -1e8 | I1_reshape(i, 2) == -1e8 | I1_reshape(i, 3) == -1e8 ...
      | I2_reshape(i, 1) == -1e8 | I2_reshape(i, 2) == -1e8 | I2_reshape(i, 3) == -1e8 )
        v(count) = i;
        count = count + 1;
    end
end

if( exist('v','var') )
    I1_reshape( v, : ) = [];
    I2_reshape( v, : ) = [];
    if use_sift(1)
        pos1(:,v) = [];
        pos2(:,v) = [];
    end
end


%% Estimate gamma values and color correction matrix H
% disp '   Estimate gamma values and matrix H ----------------------------------------------'
% if use_sift(1)
%     disp(['Number of correspondences : ',num2str(length(pos1(1,:)))]);
% end
% Initial guess for optimization step
options = optimoptions('fmincon' ,'Algorithm','sqp', 'Display', 'off');
Hl = -100.*ones(3); Hl(1,1) = 1e-4; Hl(2,2) = 1e-4; Hl(3,3) = 1e-4;
Hu =   10.*ones(3); Hu(1,1) = 1e+3; Hu(2,2) = 1e+3; Hu(3,3) = 1e+3;

min_gamma = 0.5;
max_gamma = 3.;
if isempty(ref_ref)
    H0 = eye(3);
    gamma0 = gamma_ref; 
else
    H0 = reshape(ref_ref(2:end), 3 , 3)';
    gamma0 = ref_ref(1);
end

tic;
if gamma_ref > 0
    [gammas1, fval1 ] = fmincon( @(x) myfunction11D4(x,  I1_reshape, I2_reshape  ),  [gamma0 gamma_ref H0(:)'], [], [], [], [], [min_gamma gamma_ref Hl(:)'], [max_gamma gamma_ref Hu(:)'], [], options );
    [gammas2, fval2 ] = fmincon( @(x) myfunction11D5(x,  I1_reshape, I2_reshape  ),  [gamma0 gamma_ref H0(:)'], [], [], [], [], [min_gamma gamma_ref Hl(:)'], [max_gamma gamma_ref Hu(:)'], [], options );
else 
    [gammas1, fval1 ] = fmincon( @(x) myfunction11D4(x,  I1_reshape, I2_reshape  ),  [1 1 H0(:)'], [], [], [], [], [min_gamma min_gamma Hl(:)'], [max_gamma max_gamma Hu(:)'], [], options );
    [gammas2, fval2 ] = fmincon( @(x) myfunction11D5(x,  I1_reshape, I2_reshape  ),  [1 1 H0(:)'], [], [], [], [], [min_gamma min_gamma Hl(:)'], [max_gamma max_gamma Hu(:)'], [], options );
end

%% Select the better optimization
if fval1 > fval2
    gammas = gammas2;
    H = inv(reshape(gammas2(3:end), 3 , 3)'); 
    iter = 2;
%     disp(['Error : ',num2str(fval2)]);
else
    gammas = gammas1;
    H = reshape(gammas1(3:end), 3 , 3)';
    iter = 1;
%     disp(['Error : ',num2str(fval1)]);
end
%% Check that color correction matrix is well computed
if( H(1, 1) < 0 || H(2, 2) < 0 || H(3, 3) < 0 )
    if iter>1
       H = reshape(gammas1(3:end), 3 , 3)';
    else
       H = inv(reshape(gammas2(3:end), 3 , 3)'); 
    end
    if( H(1, 1) < 0 || H(2, 2) < 0 || H(3, 3) < 0 )
        H = zeros(3,3);
    end
end
% disp('    H matrix: ' );disp(H)

% disp(['    Gamma values: ', num2str(gammas(1)), ' & ', num2str(gammas(2)) ])
% disp('    H matrix ransac: ' );disp(H) % no ransac is used for now

%% Apply gamma and H transformation to full size (original) images
%% Apply estimated gammas to the original images
l=limit^gammas(2);

% % % gammaZZZ = interp1([0 l/2 l 1],[gammas(1) gammas(1) gammas(1) 1],mean(reshape(I1, [], 3),2),'pchip');
I1(I1<0)=0;
I1_c = I1 .^ gammas(1); 

tmp0  = reshape(I1_c, [], 3);
% % % I1(I1<0) = 0;
% % % tmp0 = reshape(I1, [], 3).^gammaZZZ;

%% Apply transformation

%%

% intensity = mean(tmp0,2);
% W(intensity>=l)=0.5*(1+cos((intensity(intensity>=l) - l)*pi/(l-1)));
% W(intensity<l)=1;
% W=W.*W;
% 
% [V,D] = eig(H);
% iV = inv(V);
% 
% DW = diag(D).^W;
% % DW3x3(1,1,:)=DW(1,:);
% % DW3x3(2,2,:)=DW(2,:);
% % DW3x3(3,3,:)=DW(3,:);
% DW3(:,1,:) = DW;
% DW3(:,2,:) = DW;
% DW3(:,3,:) = DW;
% Vn = repmat(V(:),[1 length(intensity)]);
% DWiVn = DW3.*repmat(iV,[1 1 length(intensity)]);
% 
% tr = trace(H);
% 
% HW = real([ squeeze(sum(reshape(Vn.*(kron(squeeze(DWiVn(:,1,:))',ones(1,3))'),[3 3 length(intensity)]),2))' ...
%             squeeze(sum(reshape(Vn.*(kron(squeeze(DWiVn(:,2,:))',ones(1,3))'),[3 3 length(intensity)]),2))' ...
%             squeeze(sum(reshape(Vn.*(kron(squeeze(DWiVn(:,3,:))',ones(1,3))'),[3 3 length(intensity)]),2))' ]');
% HW = tr./(HW(1,:)+HW(5,:)+HW(9,:)).*HW;
% % I12tmp = W'.*I12tmp + (1-W').*tmp0;
% 
% I12 = reshape( squeeze(sum(reshape(HW.*([tmp0 tmp0 tmp0]'),[3 3 length(intensity)]),2))' , size(I1) );
% I12(I12 < 0)    = 0;
% % I12(I12 > 1)    = 1;
%%
% I12tmp = ( H * tmp0' )' ;
% intensity = mean(tmp0,2);
% W(intensity>=l)=0.5*(1+cos((intensity(intensity>=l) - l)*pi/(l-1)));
% W(intensity<l)=1;
% W=W'.*W';
% W = interp1([0 l/2 l 1],[1 1 1 0],intensity,'spline');
% I12tmp = W.*I12tmp + mean(max(I12tmp))/mean(max(tmp0))*(1-W).*tmp0;
% I12 = reshape( I12tmp, size(I1) );
% I12(I12 < 0)    = 0;
%% Interpolation of the Homography

% Id = diag( [ max(max(I2tmp(:,:,1)))/max(max(I1tmp(:,:,1))) ...
%              max(max(I2tmp(:,:,2)))/max(max(I1tmp(:,:,2))) ...
%              max(max(I2tmp(:,:,3)))/max(max(I1tmp(:,:,3))) ] );
% RGB = tmp0>l;
% intensity = mean(tmp0,2).*RGB(:,1).*RGB(:,2).*RGB(:,3);
% % intensity(intensity<l)=0;
% % imshow(reshape(intensity,[1200 1920]))
% % drawnow
% 
% if l<1
%     HW=zeros(9,length(intensity));
%     for k = 1:9
%         HW(k,:)=interp1([0 l/2 l 1],[H(k) H(k) H(k) Id(k)],intensity,'pchip');
%     end
%     I12 = reshape( squeeze(sum(reshape(HW.*(kron(tmp0,ones(1,3))'),[3 3 length(intensity)]),2))' , size(I1) );
% else
%     I12 = reshape( ( H * tmp0' )', size(I1) );
% end
% 
% I12(I12 < 0)    = 0;
% I12(I12 > 1)    = 1;

%% Interpolation of saturated pixels to bring them back in the color space

% I12 = ( H * tmp0' )' ;
% 
% saturind = I12 > 1;
% saturind = saturind(:,1)|saturind(:,2)|saturind(:,3);
% HoneLAB = rgb2lab((H * ones(3,1))');
% I12LAB = rgb2lab(I12);
% I12LAB(saturind,1) = interp1([min(unique(I12LAB(~saturind,1)')) max(unique(I12LAB(~saturind,1)')) HoneLAB(1)],[min(unique(I12LAB(~saturind,1)')) max(unique(I12LAB(~saturind,1)')) 100],I12LAB(saturind,1)');
% I12LAB(saturind,2) = interp1([min(unique(I12LAB(~saturind,2)')) max(unique(I12LAB(~saturind,2)')) HoneLAB(2)],[min(unique(I12LAB(~saturind,2)')) max(unique(I12LAB(~saturind,2)'))   0],I12LAB(saturind,2)');
% I12LAB(saturind,3) = interp1([min(unique(I12LAB(~saturind,3)')) max(unique(I12LAB(~saturind,3)')) HoneLAB(3)],[min(unique(I12LAB(~saturind,3)')) max(unique(I12LAB(~saturind,3)'))   0],I12LAB(saturind,3)');
% I12 = lab2rgb(I12LAB);
% figure
% imshow(reshape(I12,size(I1)))

%% Detection of the OER

% I12tmp = ( H * tmp0' )' ;
% I12tmp = (1-OER).*I12tmp + mean(max(I12tmp))/mean(max(tmp0))*OER.*tmp0;

% Weighted identity POSSIBLE ISSUE : Maybe I shouldn't normalize by the max
% of each channel, but rather by RGB values of the pixel with the highest intensity

Id = diag( [ max(max(I2tmp(:,:,1)))/max(max(I1tmp(:,:,1))) ...
             max(max(I2tmp(:,:,2)))/max(max(I1tmp(:,:,2))) ...
             max(max(I2tmp(:,:,3)))/max(max(I1tmp(:,:,3))) ] );
         
         
if l<1
    tmpLAB = rgb2lab(tmp0);
    
% %     Investigation in progress to improve to correction of the blacks

%     negtmpLAB = rgb2lab(1-tmp0);
%     negOER = 0.5*(tanh(1/60.*((negtmpLAB(:,1)-80)+(40-sqrt(negtmpLAB(:,2).^2+negtmpLAB(:,3).^2))))+1);
% %    0.85 seems to be a nice value 
%     negOER((negOER-0.85)<=0)=0;
%     imagesc(reshape(negOER,[1200 1920]))
    
    
    % Formula to detect the over exposed regions, 1/60,80 and 40 can be
    % customized to make the correction less "agressive"
    OER = 0.5*(tanh(1/60.*((tmpLAB(:,1)-80)+(40-sqrt(tmpLAB(:,2).^2+tmpLAB(:,3).^2))))+1);
    OER((OER-l)<=0)=0;
    
    % OER of rgb2lab([1 1 1])
    Moer = 0.880797077977882;
    
    % Normalization
    OER = 1/Moer*OER;
    
    HW=zeros(9,length(OER));
    for k = 1:9
        % Interpolation
        HW(k,:)=interp1([0 l/2/Moer l/Moer 1],[H(k) H(k) H(k) Id(k)],OER,'pchip'); 
    end
    % Applying the customized transformation with a little trick to make it
    % way quicker than a for loop
    I12 = reshape( squeeze(sum(reshape(HW.*(kron(tmp0,ones(1,3))'),[3 3 length(OER)]),2))' , size(I1) ); 
else
    I12 = reshape( ( H * tmp0' )', size(I1) );
end
% Clipping
I12(I12 < 0)    = 0;
I12(I12 > 1)    = 1;
%% Define mask
% [r_l1tmp, r_u1tmp] = discard_saturated( reshape(I1, [], 3), clip_v );
% 
% I12tmp = ( H * tmp0' )';
% I12tmp_l = I12tmp;
% I12tmp_u = I12tmp;
% for i=1:length(r_l1tmp)
%     if r_l1tmp(i)==1
%         I12tmp_l( i, : ) = -100.*ones(1, 3); 
%         I12tmp( i, : )   = -100.*ones(1, 3); 
%     end
% end
% 
% for i=1:length(r_u1tmp)
%     if r_u1tmp(i)==1
%         I12tmp_u( i, : ) = -100.*ones(1, 3); 
%         I12tmp( i, : )   = -100.*ones(1, 3); 
%     end
% end
% I12tmp = reshape( ( I12tmp' )', size(I1));
% 
% mask  = (I12tmp);
% mask( I12tmp > 0 )  = 1;
% mask  = double( mask );
% 
% I12tmp_l = reshape( ( I12tmp_l' )', size(I1));
% mask_l  = (I12tmp_l);
% mask_l( I12tmp_l > 0 )  = 1;
% mask_l  = double( mask_l );
% 
% I12tmp_u = reshape( ( I12tmp_u' )', size(I1));
% mask_u  = (I12tmp_u);
% mask_u( I12tmp_u > 0 )  = 1;
% mask_u  = double( mask_u );
% 
% mask(mask < 0) = 0;
% mask_l(mask_l < 0) = 0;
% mask_u(mask_u < 0) = 0;

%% Store in a structure gamma values and the transformation H
values.I1exp0   = I12;
% values.mask     = mask;
% values.mask_l   = mask_l;
% values.mask_u   = mask_u;
values.H        = H;
values.gamma    = gammas;
values.max_clip = ( H * ( (clip_v(2).^gammas(1)) .* ones(3,1) ) );
values.min_clip = ( H * ( (clip_v(1).^gammas(1)) .* ones(3,1) ) );


% pts1 = I1_reshape;
% pts2 = I2_reshape;
% figure;
% subplot 221
% for i =1:size(pts1, 1)
%     scatter3(pts1(i, 1), pts1(i, 2), pts1(i, 3), [], [pts1(i, 1), pts1(i, 2), pts1(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])
% subplot 222
% for i =1:length(pts2)
%     scatter3(pts2(i, 1), pts2(i, 2), pts2(i, 3), [], [pts2(i, 1), pts2(i, 2), pts2(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])
% 
% I1_c = pts1 .^ gammas(1);
%% Set to 0 negative values
%% Apply transformation
% pts12 = ( H * I1_c' )';
% pts12(pts12 < 0)   = 0;
% 
% I2_c = pts2 .^ gammas(2);
%% Set to 0 negative values
%% Apply transformation
% pts22 = ( H \ I2_c' )';
% pts22(pts22 < 0)   = 0;
% 
% pts12_c = pts12.^(1/gammas(2));
% pts12_c( pts12_c < 0 ) = 0;
% pts12_c( pts12_c > 1 ) = 1;
% subplot 224
% for i =1:size(pts12, 1)
%     scatter3(pts12(i, 1), pts12(i, 2), pts12(i, 3), [], [pts12_c(i, 1), pts12_c(i, 2), pts12_c(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])

% pts22_c = pts22.^(1/gammas(1));
% pts22_c( pts22_c < 0 ) = 0;
% pts22_c( pts22_c > 1 ) = 1;
% subplot 223
% for i =1:length(pts22)
%     scatter3(pts22(i, 1), pts22(i, 2), pts22(i, 3), [], [pts22_c(i, 1), pts22_c(i, 2), pts22_c(i, 3)], 'filled', 'MarkerEdgeColor', 'k');hold on
% end
% axis([0 1 0 1 0 1])

clearvars -except values
