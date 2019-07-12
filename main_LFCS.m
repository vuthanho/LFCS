%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Color Stabilization for Lightfield data Algorithm 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = strcat(pwd,'/../LFColorSample/take4_5');
subdir = dir(folder);
saved_data = '/saved_data/';
output = '/output/nogopt/';
%% save outputs in a .mat file
if isempty(dir(strcat(pwd,saved_data)))
    mkdir(pwd,saved_data);
end
save_file = strcat(pwd,saved_data);

%% save images as .exr files
if isempty(dir(strcat(pwd,output)))
    mkdir(pwd,output);
end
save_im = strcat(pwd,output);

%% Set up parameters for the algorithm
options.save_file = []; % name of the folder to save results
options.save_im   = save_im; % name of the folder to save images
options.exp0      = -1; % if exp0 < 0, the algorithm takes the middle exposure at the center
options.clipping  = [15 240]; % min & max values for clipping e.g. [0 255]
options.factor    = [1 0.5]; % first value is for final result, second value for computational purposes.
options.spread    = 1; % enables spreading references, 0 single ref at the center
options.sift      = 1; 
options.dense     = 0;
options.id_lut    = [];%strcat(pwd,'/src/lut/id_lut.mat');
options.limit     = 0.7; % limit after which the homography will be less applied (in linear color space)
      
startTime = tic;
lfcs_method( strcat(folder, '/'), options );
% prefixe = strcat(pwd,'/output/spreadvsnospread/0');
% for i=0.6:0.1:1
%     options.factor(2) = sqrt(i);
%     % spread
%     options.spread = 1;
%     options.save_im = strcat(prefixe,num2str(10*i),'0/spread/');
%     if isempty(dir(options.save_im))
%         mkdir(options.save_im);
%     end
%     lfcs_method( strcat(folder, '/'), options );
%     % nospread
%     options.spread = 0;
%     options.save_im = strcat(prefixe,num2str(10*i),'0/nospread/');
%     if isempty(dir(options.save_im))
%         mkdir(options.save_im);
%     end
%     lfcs_method( strcat(folder, '/'), options );
% end
t = toc( startTime );
disp(['   It took ' , num2str(t), ' s'])