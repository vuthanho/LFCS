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
names = natsort({dcm(3:end).name});
images = cell(1,64);
for i = 1:size(names,2)
        
    filename = strcat(folder,names{i});
    if (strcmp(filename(end-2:end),'exr'))
        im = exrread(filename);
        im = uint8(im.*255);
    else
        im = imread(filename);
    end
    
    images{i} = im;

end
