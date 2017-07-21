% This is a demo code for our OCT denoising algorithm:

%==========================================================================
% Read Data
%==========================================================================
clear; close all
%%% choose which example sample to use
% imgPath = '../Data/phantom_noisy.tif';
% imgPath = '../Data/orange_noisy.tif';
% imgPath = '../Data/chickenskin_noisy.tif';
imgPath = '../Data/biofilm_noisy.tif';


%%% the input image intensity should be within [0, 1]
D_input = imread(imgPath);
D_input = im2double(D_input);

% create a colormap
map = parula(1000); 
map = map(1:999,:);
map = [map;[1,1,1]];

%==========================================================================
% Statistical parameter estimation
%==========================================================================

alpha = estimatePar(D_input);

%==========================================================================
% Run Proposed Method
%==========================================================================
if ~exist('alpha','var')
    % default value
    alpha = 0.525;
end
% compute distribution coefficient
c1 = (1-alpha^2/2)^(1/4);
c2 = 1-(1-alpha^2/2)^(1/2);

% parameter selection
par.lambda = 0.4;
par.gamma = 2;
par.delta = 0.002;
par.theta = 0.98;
par.c1 = c1;
par.c2 = c2;
par.maxIter = 20;

tic;
[ U_ours_huberTV ] = ladexp_huberTV( D_input, par );
toc;

%==========================================================================
% Display images
%==========================================================================
figure; imshow(D_input,[],'Colormap',map); colorbar;
title('Noisy');
figure; imshow(U_ours_huberTV,[],'Colormap',map); colorbar;
title('Denoised');

