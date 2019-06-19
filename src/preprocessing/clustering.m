function [newdata, outliers] = clustering( data, zFactor, factor )


% data = the variable name of your array
stdDev = std(data(:)); % Compute standard deviation
meanValue = median(data(:)); % Compute mean

% Create a binary map of where outliers live.
outliers = abs(data-meanValue) > (zFactor * stdDev); 

newdata = data;
for i=1:size(outliers, 1)
    if( outliers(i, 1) == 1 | outliers(i, 2) == 1 | outliers(i, 3) == 1 )
        newdata(i, :) = factor.*ones(1, size(data, 2));
    end
end

% figure;scatter(0:size(data, 1)-1, data, 10, [1. 0 0]);hold on
% scatter(0:size(newdata, 1)-1, newdata, 10, [0. 1. 0]);