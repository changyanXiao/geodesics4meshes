clear all; close all; clc;
%--------------------------------------------------------------------------
library_directory = '../mex';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
%%
% load mask
P = imread('sala.png');
n = size(P,1);
m = size(P,2);

SE=strel('square',3);
P = imopen(P,SE);    % IMPORTANT : le bord du domaine doit pouvoir s'extraire en 4-conexité

% pick points
S = zeros(n,m);
S = int32(S);

% apply mask: IMPORTANT FOR COMPUTATIONS
S(P==0) = -1;       % negative values in seed matrix are considered outside the domain

% tensor field
T = zeros(n,m,2,2);
T(:,:,1,1) = 1;
T(:,:,2,2) = 1;

% sampling
nb_seeds = 100;
tic
[U,V,dUx,dUy] = anisoFPMeshing2Diterative(T,S,nb_seeds);   % 5th USE (with mask implicitely in S)
toc
%%
% display Voronoi diagram
%clf(ctrl);
figure;
U(U==-1)=-10;
sp=[];
for i=1:n
    for j=1:m
        if (S(i,j)>0)
            sp=[sp,[i;j]];
        end
    end
end
imagesc(U);
axis image; axis off;
colormap jet(256);
hold on;
plot_pts(sp, 'r.',10);
hold off;
figure;
imagesc(V);
axis image; axis off;
hold on;
plot_pts(sp,'k.',20);
hold off;
