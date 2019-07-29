function [] = lfcs_method( folder, options )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculate color stabilized images of the lightfield array
%% Set a reference and estimate gamma and color
%%
%% Inputs:  1. folder  -> char name of the folder where LF array is
%%          2. options -> struct containing save_file/save_im/exp0/clipping
%%                          /factor/sift/dense/spread/id_lut/limit/
%%
%% Outputs: 1. None
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check input values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2
    %% Set up parameters for the algorithm
    
    if isa(options,'char') % If the algorithm is compiled as an executable
        indpar = find(options==',');
        options = strcat(options(1),'''',options(2:(indpar(1)-1)),'''',... % save_file
                    options(indpar(1)),'''',options((indpar(1)+1):(indpar(2)-1)),'''',... % save_im
                    options(indpar(2):end)); % Yeah this is ugly but necessary and still needs to be done for options{9} too (location of the LUT)
        options = eval(options);
        save_file = options{1};
        save_im   = options{2};
        exp0      = options{3};
        clipping  = options{4};
        factor    = options{5};
        use_sift  = [options{6} options{7}];
        write     = ~isempty(options{2});
        spread    = options{8};
        if ~isempty(options{9})
            id_lut    = cell2mat(struct2cell(load(options{9})));
        else
            id_lut    = idlutmaker(33);
        end
        l     = options{10};
        if isempty(save_file)
            to_save = 0;
        else
            to_save = 1;
        end
    else 
        save_file = options.save_file;
        save_im   = options.save_im;
        exp0      = options.exp0;
        clipping  = options.clipping;
        factor    = options.factor;
        use_sift  = [options.sift options.dense];
        write     = ~isempty(options.save_im);
        spread    = options.spread;
        if ~isempty(options.id_lut)
            id_lut    = cell2mat(struct2cell(load(options.id_lut)));
        else
            id_lut    = idlutmaker(33);
        end
        l     = options.limit;
        if isempty(save_file)
            to_save = 0;
        else
            to_save = 1;
        end
    end
    
    
end
if nargin < 2
    
    %% Set up parameters for the algorithm
    if ismac
        save_file = strcat(pwd,'/CUBEs/');
    elseif isunix
        save_file = strcat(pwd,'/CUBEs/');
    elseif ispc
        save_file = strcat(pwd,'\CUBEs\');
    else
        disp('Platform not supported')
    end
    save_im   = [];
    exp0      = -1;
    to_save   = 1;
    clipping  = [15 240];
    factor    = [1 .5];
    use_sift  = [1 0];
    write     = 0;
    spread    = 1;
    id_lut    = idlutmaker(33);
    l     = 0.75;
end
if ~isempty(save_file)
    if isempty(dir(save_file))
        mkdir(save_file);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read images and exposures
disp 'Reading images ...'
im = read_images( folder );
% imnc = read_images( '/home/bonnetvalentin/Documents/LFCS/LFColorSample/take4_5/' ); %% Only to apply the results of the cropped color checker to non cropped views

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
    values{i}.H        = zeros(3, 3);
    values{i}.gamma    = zeros(1, 2);
    values{i}.max_clip = zeros(1, 3);
    values{i}.min_clip = zeros(1, 3);
end

%% Initialize reference exposure image

im_med_exp0 = ( double( imresize(im{exp0}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );


%% gamma selection

gamma = 1;

%% First calculation for the center

for i = center
    disp(['Calculating image ', num2str(i),' --> image ', num2str(exp0)  ])

    im_med = ( double( imresize(im{i}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );

%     Calculating the homography
    values{i} = compute_colorstabilization( double( imresize(im{i}, factor(1),'bilinear', 'Colormap', 'original') )./max_im, ...
        im_med, im_med_exp0, clip, gamma, use_sift, [] , l);
    values{i}.I1exp0 = (values{i}.I1exp0).^(1/values{i}.gamma(2));
%     Saving the LUT
    if to_save
        lut = id_lut .^ values{i}.gamma(1);
        H = values{i}.H;
        limit=l^values{i}.gamma(2);

% %         Detection of overexposured region

        Id = diag( [ max(max(im_med_exp0(:,:,1)))/max(lut(:,1)) max(max(im_med_exp0(:,:,2)))/max(lut(:,2)) max(max(im_med_exp0(:,:,3)))/max(lut(:,3)) ] );
        if limit<1
            tmpLAB = rgb2lab(lut);
            OER = 0.5*(tanh(1/60.*((tmpLAB(:,1)-80)+(40-sqrt(tmpLAB(:,2).^2+tmpLAB(:,3).^2))))+1);
            OER((OER-limit)<=0)=0;
            Moer = 0.880797077977882; % OER of rgb2lab([1 1 1])
            OER = 1/Moer*OER; 
            HW=zeros(9,length(OER));
            for k = 1:9
                HW(k,:)=interp1([0 limit/2/Moer limit/Moer 1],[H(k) H(k) H(k) Id(k)],OER,'pchip');
            end
            lut = reshape( squeeze(sum(reshape(HW.*(kron(lut,ones(1,3))'),[3 3 length(OER)]),2))' , size(lut) );
        else
            lut = ( H * lut' )';
        end
        lut(lut < 0)    = 0;
        lut(lut > 1)    = 1;
        lut = lut.^(1/values{i}.gamma(2));
        write_cube(strcat(save_file,num2str(i),'.CUBE'),num2str(i), [0.0 0.0 0.0], [1.0 1.0 1.0], lut' );
%         struct_name = strcat(save_file,'lut',num2str(i), '.mat');
%         save(struct_name, 'lut', '-v7.3');
    end 

% %     Only to apply the results of the cropped color checker to non
%       cropped views
%     H = values{i}.H;
%     limit=l^values{i}.gamma(2);
%     lut = reshape(im2double(imnc{i}), [], 3);
%     Id = diag( [ max(max(im_med_exp0(:,:,1)))/max(lut(:,1)) max(max(im_med_exp0(:,:,2)))/max(lut(:,2)) max(max(im_med_exp0(:,:,3)))/max(lut(:,3)) ] );
%     if limit<1
%         tmpLAB = rgb2lab(lut);
%         OER = 0.5*(tanh(1/60.*((tmpLAB(:,1)-80)+(40-sqrt(tmpLAB(:,2).^2+tmpLAB(:,3).^2))))+1);
%         OER((OER-limit)<=0)=0;
%         Moer = 0.880797077977882; % OER of rgb2lab([1 1 1])
%         OER = 1/Moer*OER; 
%         HW=zeros(9,length(OER));
%         for k = 1:9
%             HW(k,:)=interp1([0 limit/2/Moer limit/Moer 1],[H(k) H(k) H(k) Id(k)],OER,'pchip');
%         end
%         lut = reshape( squeeze(sum(reshape(HW.*(kron(lut,ones(1,3))'),[3 3 length(OER)]),2))' , size(lut) );
%     else
%         lut = ( H * lut' )';
%     end
%     lut(lut < 0)    = 0;
%     lut(lut > 1)    = 1;
%     lut = lut.^(1/values{i}.gamma(2));
%     lut = reshape(lut,size(imnc{i}));
%     if spread
%         filename = strcat('/home/bonnetvalentin/Documents/LFCS/Test/output/cropped',num2str(round(10*(factor(2)^2))),'0/',num2str(i),'.exr');
%     else
%         filename = strcat('/home/bonnetvalentin/Documents/LFCS/Test/output/center_cropped',num2str(round(10*(factor(2)^2))),'0/',num2str(i),'.exr');
%     end
%     exrwrite(lut,filename);

    
    %% Define the limit ranges for the weighting function
    clip_vR(i,:) = [ values{i}.min_clip(1) values{i}.max_clip(1)];
    clip_vG(i,:) = [ values{i}.min_clip(2) values{i}.max_clip(2)];
    clip_vB(i,:) = [ values{i}.min_clip(3) values{i}.max_clip(3)];
    
    %% Display the range values
%     disp(['Clipping min [R, G, B]: [',num2str( clip_vR(i,1) ), ', ' num2str( clip_vG(i,1) ), ', ', num2str( clip_vB(i,1) ), ']' ])
%     disp(['Clipping max [R, G, B]: [',num2str( clip_vR(i,2) ), ', ' num2str( clip_vG(i,2) ), ', ', num2str( clip_vB(i,2) ), ']' ])
    
    clear im_med
end

%% Number of contours from the center
if isPeven
    nbrContours = p/2-1;
else
    nbrContours = floor(p/2);
end

%% Whole calculation

for nbr = 1:nbrContours
    if isPeven
        Contour = getContourE(p,nbr);
    else
        Contour = getContourO(p,nbr);
    end
    for image = Contour
        i = image(1);
%         exp_neighbour = image(2);
        
        im_med = ( double( imresize(im{i}, factor(2),'bilinear', 'Colormap', 'original') )./max_im );

        %% Calculation
        if ~spread
            disp(['Calculating image ', num2str(i),' --> image ', num2str(exp0)  ])
            values{i} = compute_colorstabilization( double( imresize(im{i},factor(1),'bilinear' , 'Colormap', 'original') )./max_im, ...
            im_med, im_med_exp0,... %imresize(values{exp0}.I1exp0, factor(2),'bilinear', 'Colormap', 'original'),  ...
            clip, gamma, use_sift,[],l);%[values{exp_neighbour}.gamma(1) values{exp_neighbour}.H(:)'] );
        else
            disp(['Calculating image ', num2str(i),' --> image ', num2str(image(2))])
            values{i} = compute_colorstabilization( double( imresize(im{i},factor(1),'bilinear' , 'Colormap', 'original') )./max_im, ...
            im_med, imresize(values{image(2)}.I1exp0, factor(2),'bilinear', 'Colormap', 'original'),  ...
            clip, gamma, use_sift,[],l);%[values{exp_neighbour}.gamma(1) values{exp_neighbour}.H(:)'] );
        end
        values{i}.I1exp0 = (values{i}.I1exp0).^(1/values{i}.gamma(2));
%     Saving the LUT
        if to_save
            lut = id_lut .^ values{i}.gamma(1);
            H = values{i}.H;        
            limit=l^values{i}.gamma(2);

% %             Detection of OER

            Id = diag( [ max(max(values{image(2)}.I1exp0(:,:,1)))/max(lut(:,1)) ...
                         max(max(values{image(2)}.I1exp0(:,:,2)))/max(lut(:,2)) ...
                         max(max(values{image(2)}.I1exp0(:,:,3)))/max(lut(:,3)) ] );
            if limit<1
                tmpLAB = rgb2lab(lut);
                OER = 0.5*(tanh(1/60.*((tmpLAB(:,1)-80)+(40-sqrt(tmpLAB(:,2).^2+tmpLAB(:,3).^2))))+1);
                OER((OER-limit)<=0)=0;
                Moer = 0.880797077977882; % OER of rgb2lab([1 1 1])
                OER = 1/Moer*OER; 
                HW=zeros(9,length(OER));
                for k = 1:9
                    HW(k,:)=interp1([0 limit/2/Moer limit/Moer 1],[H(k) H(k) H(k) Id(k)],OER,'pchip');
                end
                lut = reshape( squeeze(sum(reshape(HW.*(kron(lut,ones(1,3))'),[3 3 length(OER)]),2))' , size(lut) );
            else
                lut = ( H * lut' )';
            end
            lut(lut < 0)    = 0;
            lut(lut > 1)    = 1;
            lut = lut.^(1/values{i}.gamma(2));
            write_cube(strcat(save_file,num2str(i),'.CUBE'),num2str(i), [0.0 0.0 0.0], [1.0 1.0 1.0], lut' );
%             struct_name = strcat(save_file,'lut',num2str(i), '.mat');
%             save(struct_name, 'lut', '-v7.3');
        end
        
% %     Only to apply the results of the cropped color checker to non
%       cropped views
%         H = values{i}.H;
%         limit=l^values{i}.gamma(2);
%         lut = reshape(im2double(imnc{i}), [], 3);
%         Id = diag( [ max(max(im_med_exp0(:,:,1)))/max(lut(:,1)) max(max(im_med_exp0(:,:,2)))/max(lut(:,2)) max(max(im_med_exp0(:,:,3)))/max(lut(:,3)) ] );
%         if limit<1
%             tmpLAB = rgb2lab(lut);
%             OER = 0.5*(tanh(1/60.*((tmpLAB(:,1)-80)+(40-sqrt(tmpLAB(:,2).^2+tmpLAB(:,3).^2))))+1);
%             OER((OER-limit)<=0)=0;
%             Moer = 0.880797077977882; % OER of rgb2lab([1 1 1])
%             OER = 1/Moer*OER; 
%             HW=zeros(9,length(OER));
%             for k = 1:9
%                 HW(k,:)=interp1([0 limit/2/Moer limit/Moer 1],[H(k) H(k) H(k) Id(k)],OER,'pchip');
%             end
%             lut = reshape( squeeze(sum(reshape(HW.*(kron(lut,ones(1,3))'),[3 3 length(OER)]),2))' , size(lut) );
%         else
%             lut = ( H * lut' )';
%         end
%         lut(lut < 0)    = 0;
%         lut(lut > 1)    = 1;
%         lut = lut.^(1/values{i}.gamma(2));
%         lut = reshape(lut,size(imnc{i}));
%         if spread
%             filename = strcat('/home/bonnetvalentin/Documents/LFCS/Test/output/cropped',num2str(round(10*(factor(2)^2))),'0/',num2str(i),'.exr');
%         else
%             filename = strcat('/home/bonnetvalentin/Documents/LFCS/Test/output/center_cropped',num2str(round(10*(factor(2)^2))),'0/',num2str(i),'.exr');
%         end
%         exrwrite(lut,filename);
        
        %% Define the limit ranges for the weighting function
        clip_vR(i,:) = [ values{i}.min_clip(1) values{i}.max_clip(1)];
        clip_vG(i,:) = [ values{i}.min_clip(2) values{i}.max_clip(2)];
        clip_vB(i,:) = [ values{i}.min_clip(3) values{i}.max_clip(3)];

        %% Display the range values
%         disp(['Clipping min [R, G, B]: [',num2str( clip_vR(i,1) ), ', ' num2str( clip_vG(i,1) ), ', ', num2str( clip_vB(i,1) ), ']' ])
%         disp(['Clipping max [R, G, B]: [',num2str( clip_vR(i,2) ), ', ' num2str( clip_vG(i,2) ), ', ', num2str( clip_vB(i,2) ), ']' ])

        clear im_med
    end
end

%% Write the images to the specified directory if it's enabled

if write
    if isempty(dir(save_im))
        mkdir(save_im);
    end
    for i = 1:size(values,2)
    filename = strcat(save_im,num2str(i),'.exr');
    exrwrite(values{i}.I1exp0,filename);
    end
end

clear all