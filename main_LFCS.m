%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Color Stabilization for Lightfield data Algorithm 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = strcat(pwd,'/../LFColorSample/take2_1');
subdir = dir(folder);
saved_data = '/saved_data/';
output = '/output/PNG/take2_1/';
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

%% Set up the parameters for RANSAC
coef.minPtNum    = 5;
coef.iterNum     = 200;
coef.thDist      = 0.005;  
coef.thInlrRatio = 0.05; 

%% Set up parameters for the algorithm
options.save_file = save_file; % name of the folder to save results
options.save_im   = save_im; % name of the folder to save images
options.write     = 1; % 1 enables save image 
options.show      = 0; % 1 show / 0 not show
options.exp0      = 29; % if exp0 < 0, the algorithm takes the middle exposure at the center
options.clipping  = [15 240]; % min & max values for clipping e.g. [0 255]
options.factor    = [1 .5]; % first value is for final result, second value for computational purposes.
options.spread    = 1; % enables spreading references, 0 single ref at the center
options.sift      = 1; 
options.prog      = 0;
options.dense     = 0;
options.homography = 0; % 0 3x3 / 1 affine / 2 homography
options.id_lut    = strcat(pwd,'/src/lut/id_lut33.mat');
      
startTime = tic;

lfcs_method( strcat(folder, '/'), coef, options );

t = toc( startTime );
disp(['   It took ' , num2str(t), ' s'])