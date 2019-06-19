function lab_pts = rgb2lab_pts( pts )

rgb2xyz = [0.4124564  0.3575761  0.1804375;0.2126729  0.7151522  0.0721750;0.0193339  0.1191920  0.9503041];

xyz_pts = (rgb2xyz * pts')';

white_point = [95.047,100.000,108.883]./100;

lab_pts = changeXYZ2LAB(xyz_pts, white_point);

lab_pts(:, 1) =  lab_pts(:, 1)./100;
% lab_pts(:, 2) = (lab_pts(:, 2) + 128) ./ 255;
% lab_pts(:, 3) = (lab_pts(:, 3) + 128) ./ 255;
lab_pts(:, 2) = (lab_pts(:, 2) +   86.185) ./ 184.439;
lab_pts(:, 3) = (lab_pts(:, 3) +  107.863) ./ 202.345;
