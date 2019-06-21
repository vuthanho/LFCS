function [] = lfcs_method( folder, coef, options )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculate color stabilized images of the lightfield array
%% Set a reference and estimate gamma and color
%%
%% Inputs:  1. folder  -> char name of the folder where LF array is
%%          2. coef    -> struct for RANSAC parameters, (not used for now)
%%          3. options -> struct containing save_file/show/exp0/clipping
%%                               factor/sift/dense
%%
%% Outputs: 1. None
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check input values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 3
    
    %% Set up parameters for the algorithm
    save_file = options.save_file;
    save_im   = options.save_im;
    show      = options.show;
    exp0      = options.exp0;
    clipping  = options.clipping;
    factor    = options.factor;
    use_sift  = [options.sift options.dense];
    hom       = options.homography;
    write     = options.write;
    spread    = options.spread;
    prog      = options.prog;
    id_lut    = cell2mat(struct2cell(load(options.id_lut)));
    limit     = options.limit;
    if isempty(save_file)
        to_save = 0;
    else
        to_save = 1;
    end
    
end
if nargin < 3
    
    %% Set up parameters for the algorithm
    save_file = [];
    to_save   = 0;
    show      = 0;
    clipping  = [0 255];
    factor    = [1 .25];
    use_sift  = [0 0];
end
if nargin < 2
    
    %% Set up the parameters for RANSAC
    coef.minPtNum    = 5;
    coef.iterNum     = 100;
    coef.thDist      = 0.09;
    coef.thInlrRatio = .001;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read images and exposures
disp 'Reading images ...'
im = read_images( folder );

P = length(im); % number of photos
p = sqrt(P);
isPeven = (mod(p,2)==0);
center = getCenter(p);
%% Select the reference exposure among the center images if needed
if exp0 < 0
    if isPeven
        exp0 = center(choose_best_exposure( im(center), clipping )); % best among the four at the center
    else
        exp0 = center(1); % just select the center
    end
end

%% Initialize clip values vector
clip_vR = zeros( P, 2);
clip_vG = zeros( P, 2);
clip_vB = zeros( P, 2);

clip = clipping ./ 255;

%% Check image class
if isa( im{1}, 'uint8' )
    max_im = 255.;
end
if isa( im{1}, 'uint16' )
    max_im = 65535.;
end

%% Allocate memory for outputs
values = cell(1,P);
output_size = size(imresize(im{exp0},factor(1),'bilinear') );
for i=1:P
    values{i}.I1exp0   = zeros(output_size);
%     values{i}.mask     = zeros( size(im{exp0}, 1), size(im{exp0}, 2) );
    values{i}.H        = zeros(3, 3);
    values{i}.gamma    = zeros(1, 2);
    values{i}.max_clip = zeros(1, 3);
    values{i}.min_clip = zeros(1, 3);
end

%% Initialize reference exposure image
if prog
    im_med_exp0 = ( double( imresize(im{exp0}, prog,'bilinear', 'Colormap', 'original') )./max_im );
else
    im_med_exp0 = ( double( imresize(im{exp0}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );
end

%% Best gamma selection
med_value = median( im_med_exp0(:) );
if( med_value < 0.33 )
    gamma = 1.8;
elseif( med_value >= 0.33 && med_value < 0.66)
    gamma = 2.;
else
    gamma = 2.2;
end
gamma = 1;
% gamma = 1; %% Test with forced gamma = 1 and without
disp([' Gamma value ', num2str(gamma)])


%% First calculation for the center

for i = center
    disp(['Calculating image ', num2str(i),' --> image ', num2str(exp0)  ])
    if prog
        im_med = ( double( imresize(im{i}, prog,'bilinear', 'Colormap', 'original') )./max_im );
    else
        im_med = ( double( imresize(im{i}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );
    end
% % % %     gamma = -1;
    switch( hom )
        case 0
            values{i} = compute_colorstabilization( double( imresize(im{i}, factor(1),'bilinear', 'Colormap', 'original') )./max_im, ...
                im_med, im_med_exp0,  ...
                coef, clip, gamma, use_sift, [] , limit);
        case 1
            values{i} = compute_colorstabilization_aff( double( imresize(im{i}, factor(1), 'Colormap', 'original') )./max_im, ...
                im_med, im_med_exp0,  ...
                coef, clip, gamma, use_sift );
        case 2
            values{i} = compute_colorstabilization_hom( double( imresize(im{i}, factor(1), 'Colormap', 'original') )./max_im, ...
                im_med, im_med_exp0,  ...
                coef, clip, gamma, use_sift );
            
    end
    if to_save
        l=limit^values{i}.gamma(1);
        lut = id_lut .^ values{i}.gamma(1); 
        Id = mean(sum(values{i}.H,2))*eye(3);
        RGB = lut>l;
        intensity = min(1,mean(lut,2)).*RGB(:,1).*RGB(:,2).*RGB(:,3);
        HW=zeros(9,length(intensity));
        for k = 1:9
            HW(k,:)=interp1([0 l/2 l 1],[values{i}.H(k) values{i}.H(k) values{i}.H(k) Id(k)],intensity,'pchip');
        end
        lut = reshape( squeeze(sum(reshape(HW.*(kron(lut,ones(1,3))'),[3 3 length(intensity)]),2))' , size(lut) );
        lut(lut<0) = 0;
        lut(lut>1.0) = 1.0;
        write_cube(strcat(save_file,num2str(i),'.CUBE'),num2str(i), [0.0 0.0 0.0], [1.0 1.0 1.0], lut' );
        struct_name = strcat(save_file,'lut',num2str(i), '.mat');
        save(struct_name, 'lut', '-v7.3');
    end
    gamma = values{center(1)}.gamma(2);
% % % %     gamma = values{i}.gamma(2);
    values{i}.I1exp0 = (values{i}.I1exp0).^(1/gamma);
    %% Define the limit ranges for the weighting function
    clip_vR(i,:) = [ values{i}.min_clip(1) values{i}.max_clip(1)];
    clip_vG(i,:) = [ values{i}.min_clip(2) values{i}.max_clip(2)];
    clip_vB(i,:) = [ values{i}.min_clip(3) values{i}.max_clip(3)];
    
    %% Display the range values
%     disp(['Clipping min [R, G, B]: [',num2str( clip_vR(i,1) ), ', ' num2str( clip_vG(i,1) ), ', ', num2str( clip_vB(i,1) ), ']' ])
%     disp(['Clipping max [R, G, B]: [',num2str( clip_vR(i,2) ), ', ' num2str( clip_vG(i,2) ), ', ', num2str( clip_vB(i,2) ), ']' ])
    
    clear im_s im_med
end

%% Number of contours from the center
if isPeven
    nbrContours = p/2-1;
else
    nbrContours = floor(p/2);
end

%% Whole calculation

L_factor = interp1([0 nbrContours],[prog factor(2)],0:nbrContours);
for nbr = 1:nbrContours
    if prog
        factor(2) = L_factor(nbr+1);
        if ~spread
            im_med_exp0 = ( double( imresize(im{exp0}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );
        end
    end
    if isPeven
        Contour = getContourE(p,nbr);
    else
        Contour = getContourO(p,nbr);
    end
    for image = Contour
        i = image(1);
        exp_neighbour = image(2);
       
      
        %% Initialize reference exposure image
        disp(['Calculating image ', num2str(i),' --> image ', num2str(image(2)), ' | Resize = ', num2str(factor(2))])

        im_med = ( double( imresize(im{i}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );
% % % %         gamma = -1;
        %% Calculation
        switch( hom )
            case 0
                if ~spread
                    values{i} = compute_colorstabilization( double( imresize(im{i},factor(1),'bilinear' , 'Colormap', 'original') )./max_im, ...
                    im_med, im_med_exp0,... %imresize(values{exp0}.I1exp0, factor(2),'bilinear', 'Colormap', 'original'),  ...
                    coef, clip, gamma, use_sift,[]);%[values{exp_neighbour}.gamma(1) values{exp_neighbour}.H(:)'] );
                else
                    
%                     values{i} = compute_colorstabilization( double( imresize(im{i},factor(1),'bilinear' , 'Colormap', 'original') )./max_im, ...
%                         im_med, imresize(reshape((values{image(2)}.H*(reshape(double( imresize(im{image(2)},factor(1),'bilinear' , 'Colormap', 'original') )./max_im,[],3)...
%                         .^(values{image(2)}.gamma(1) / values{image(2)}.gamma(2)))')',size(values{image(2)}.I1exp0)), factor(2),'bilinear', 'Colormap', 'original'),  ...
%                         coef, clip, gamma, use_sift,[]);
                    
                    values{i} = compute_colorstabilization( double( imresize(im{i},factor(1),'bilinear' , 'Colormap', 'original') )./max_im, ...
                    im_med, imresize(values{image(2)}.I1exp0, factor(2),'bilinear', 'Colormap', 'original'),  ...
                    coef, clip, gamma, use_sift,[],limit);%[values{exp_neighbour}.gamma(1) values{exp_neighbour}.H(:)'] );
                end
            case 1
                values{i} = compute_colorstabilization_aff( double( imresize(im{i}, factor(1), 'Colormap', 'original') )./max_im, ...
                    im_med, im_med_ref,  ...
                    coef, clip, gamma, use_sift );
            case 2
                values{i} = compute_colorstabilization_hom( double( imresize(im{i}, factor(1), 'Colormap', 'original') )./max_im, ...
                    im_med, im_med_ref,  ...
                    coef, clip, gamma, use_sift );
        end
        if to_save
            l=limit^values{i}.gamma(1);
            lut = id_lut .^ values{i}.gamma(1); 
            Id = mean(sum(values{i}.H,2))*eye(3);
            RGB = lut>l;
            intensity = min(1,mean(lut,2)).*RGB(:,1).*RGB(:,2).*RGB(:,3);
            HW=zeros(9,length(intensity));
            for k = 1:9
                HW(k,:)=interp1([0 l/2 l 1],[values{i}.H(k) values{i}.H(k) values{i}.H(k) Id(k)],intensity,'pchip');
            end
            lut = reshape( squeeze(sum(reshape(HW.*(kron(lut,ones(1,3))'),[3 3 length(intensity)]),2))' , size(lut) );
            lut(lut<0) = 0;
            lut(lut>1.0) = 1.0;
            write_cube(strcat(save_file,num2str(i),'.CUBE'),num2str(i), [0.0 0.0 0.0], [1.0 1.0 1.0], lut' );
            struct_name = strcat(save_file,'lut',num2str(i), '.mat');
            save(struct_name, 'lut', '-v7.3');
        end
% % % %         gamma = values{i}.gamma(2);
        values{i}.I1exp0 = (values{i}.I1exp0).^(1/gamma);
        %% Define the limit ranges for the weighting function
        clip_vR(i,:) = [ values{i}.min_clip(1) values{i}.max_clip(1)];
        clip_vG(i,:) = [ values{i}.min_clip(2) values{i}.max_clip(2)];
        clip_vB(i,:) = [ values{i}.min_clip(3) values{i}.max_clip(3)];

        %% Display the range values
%         disp(['Clipping min [R, G, B]: [',num2str( clip_vR(i,1) ), ', ' num2str( clip_vG(i,1) ), ', ', num2str( clip_vB(i,1) ), ']' ])
%         disp(['Clipping max [R, G, B]: [',num2str( clip_vR(i,2) ), ', ' num2str( clip_vG(i,2) ), ', ', num2str( clip_vB(i,2) ), ']' ])

        clear im_s im_med
    end
end

%% Write the images to the specified directory if it's enabled

if write
    for i = 1:size(values,2)
    filename = strcat(save_im,num2str(i),'.png');
    imwrite(values{i}.I1exp0,filename);
    end
end

%% Save homographies
% if to_save
%     struct_name = strcat(save_file, 'data.mat');
%     save(struct_name, 'values', '-v7.3');
% end

clear all