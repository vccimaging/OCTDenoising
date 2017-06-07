function [ ratio ] = estimatePar( img )
%ESTIMATEPAR estimates the constant ratio between the standard deviation
%and mean of homogeneous regions in OCT data
M = prctile(img(:),99);
h = figure; imshow(img,[0 M]);

fprintf('\n...\nPlease select a homogeneous area in the figure for parameter estimation \n...\n');
% select a homogeneous area
rect = getrect(h);
region = img(rect(1): rect(1)+rect(3), rect(2): rect(2)+rect(4));

ratio = std(region(:),1)/mean(region(:));

close(h);
end

