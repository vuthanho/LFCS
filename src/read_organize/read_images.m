function images = read_images( folder )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Read images from a given folder, and store them in a cell array
%%
%%  Input:  1. folder -> char containing the name of the folder
%%
%%  Output: 1. images -> cell of P images of size mxn
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dcm = dir(folder);

for i = 3:size(dcm,1)
        
    filename = strcat(folder,dcm(i).name);
    if (strcmp(filename(end-2:end),'exr'))
        im = exrread(filename);
        im = uint8(im.*255);
    else
        im = imread(filename);
    end
    
    images{i-2} = im;

end
