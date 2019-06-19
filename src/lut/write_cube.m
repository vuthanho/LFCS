function [] = write_cube(filename, lut_title, lut_min, lut_max, lut )

% opening or writing over the file
fileID = fopen(filename,'wt');

% writting the title
fprintf(fileID,strcat('TITLE "',lut_title,'"\n'));
fprintf(fileID,'\n');

% writting the LUT size
fprintf(fileID,'LUT_%dD_SIZE %d\n',length(lut(:,1)),round(length(lut)^(1/3)));
fprintf(fileID,'\n');

% writting the domain min and max
fprintf(fileID,'DOMAIN_MIN %f %f %f\n',lut_min);
fprintf(fileID,'DOMAIN_MAX %f %f %f\n',lut_max);
fprintf(fileID,'\n');

% writting the LUT
fprintf(fileID,'%f %f %f\n',lut);

fclose(fileID);
end

