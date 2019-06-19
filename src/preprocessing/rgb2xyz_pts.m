function xyz_pts = rgb2xyz_pts( pts )

rgb2xyz = [0.4124564  0.3575761  0.1804375;0.2126729  0.7151522  0.0721750;0.0193339  0.1191920  0.9503041];

xyz_pts = (rgb2xyz * pts')';
