% This is a demo code for our denoising algorithm
clear; close all

%==========================================================================
% Read Data
%==========================================================================
%%% choose which example sample to use
% load('../Data/phantom_noisy.mat');
% load('../Data/orange_noisy.mat');
% load('../Data/chickenskin_noisy.mat');
load('../Data/biofilm_noisy.mat');

% create a colormap
map = parula(1000); 
map = map(1:999,:);
map = [map;[1,1,1]];
M = prctile(D_input(:),99);

figure; imshow(D_input,[0 M],'Colormap',map); colorbar;
title('Noisy');

%==========================================================================
% Run Proposed Method
%==========================================================================

% distribution coefficient estimation
alpha = 0.525;
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

figure;clf; imshow(U_ours_huberTV,[0 M],'Colormap',map); colorbar;
title('Denoised');


