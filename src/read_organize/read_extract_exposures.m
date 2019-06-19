function [images, exposures, name] = read_extract_exposures( folder )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Read images from a given folder, and store them in a cell array with 
%%  the info information for exposure time
%%
%%  Input:  1. folder -> char containing the name of the folder
%%
%%  Output: 1. images -> cell of P images of size mxn
%%          2. exposures -> array 1XP containing the exposure times for each image
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dcm = dir(folder);

for i = 3:size(dcm,1)
        
    filename = strcat(folder, dcm(i).name);
    im = imread(filename);
    info = imfinfo(filename);
    
    images{i-2} = im;
    exposures(i-2) = info.DigitalCamera.ExposureTime;
end