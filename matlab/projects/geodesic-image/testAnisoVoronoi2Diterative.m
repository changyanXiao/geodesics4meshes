clear all; close all; clc;
%--------------------------------------------------------------------------
library_directory = '../mex';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);

%% load image
P = imread('32.bmp');
P = double(P(:,:,1));
P = rescale(P, 0, 1);
n = size(P,1);
m = size(P,2);

%% FIRST EXAMPLE: Euclidean Voronoi diagram from 3 seeds

% compute seed matrix (with labels)
S = zeros(n,m);
S = int32(S);

ctrl = figure;
imageplot(S);
title('Pick 3 points');
axis image; axis off;
hold on;
sp = pick_points(3);
hold off;
sp = int32(sp);

label = 1;
for i=1:size(sp,2)
   S(sp(1,i),sp(2,i)) = label;
   label = label+1;
end

% tensor field
T = zeros(n,m,2,2);
T(:,:,1,1) = 1;
T(:,:,2,2) = 1;

% compute distance en Voronoi maps
tic
[U,V,dUx, dUy] = anisoVoronoi2Diterative(T,S);   % FIRST USE
toc

% display Voronoi diagram
clf(ctrl);
imagesc(V);
axis image; axis off;
hold on;
plot_pts(sp, 'r.',10);
hold off;

%% SECOND EXAMPLE: Dynamic update of the Voronoi diagram by adding seeds

figure(ctrl);
nb_seeds = 4;

for i=1:nb_seeds
   hold on;
   p =  pick_points(1);
   hold off;
   p = int32(p);
   sp = [sp,p];
   
   % important: updates needed for computations
   S(p(1),p(2)) = label;               % UPDATE OF SEEDS
   U(p(1),p(2)) = 0.0;                 % UPDATE OF DISTANCE
   V(p(1),p(2)) = label;               % UPDATE OF VORONOI
   
   label = label+1;
   
   tic
   [dUx,dUy] = anisoVoronoi2Diterative(T,S,U,V);   % SECOND USE
   toc
   
   % display Voronoi diagram
   clf(ctrl);
   imagesc(V);
   axis image; axis off;
   hold on;
   plot_pts(sp, 'r.',10);
   hold off;
end

% comparison
[U0,V0,dUx,dUy] = anisoVoronoi2Diterative(T,S);
ctrl2 = figure;
VV = V-V0;
imagesc(VV);
axis image; axis off;

%% THIRD EXAMPLE: Dynamic (local) update of the Voronoi diagram by adding seeds

figure(ctrl);
nb_seeds = 4;

for i=1:nb_seeds
   hold on;
   p =  pick_points(1);
   hold off;
   p = int32(p);
   sp = [sp,p];
   
   % important: updates needed for computations
   S(p(1),p(2)) = label;               % UPDATE OF SEEDS
   U(p(1),p(2)) = 0.0;                 % UPDATE OF DISTANCE
   V(p(1),p(2)) = label;               % UPDATE OF VORONOI
   
   label = label+1;
   
   tic
   [dUx,dUy] = anisoVoronoi2Diterative(T,S,U,V,double(p));   % THIRD USE (from p only)
   toc
   
   % display Voronoi diagram
   clf(ctrl);
   imagesc(V);
   axis image; axis off;
   hold on;
   plot_pts(sp, 'r.',10);
   hold off;
end

% comparison
[U0,V0,dUx,dUy] = anisoVoronoi2Diterative(T,S);
ctrl2 = figure;
VV = V-V0;
imagesc(VV);
axis image; axis off;
%% 4th EXAMPLE: Vornoi segmentation from scribbles

% compute structure tensor around the curve (??pas terrible trouver un bon champ??)
epsilon = 0.000001;
alpha = 1;
sigma1 = 1.4;   % for gradient computation
sigma2 = 2;    % for filtering the tensor field
options.use_anisotropic = 0;
T = compute_structure_tensor(P, sigma1, sigma2);
[e1,e2,l1,l2] = perform_tensor_decomp(T);
T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);

% get scribbles
[S,label] = interactive_labeling(P);
S = int32(S);
tic
  [U,V,dUx,dUy] = anisoVoronoi2Diterative(T,S);   % 4th USE (from contours)
toc

ctrl = clf;
M = P;
M(S>0) = mmax(P)+1;
figure(ctrl);
imagesc(M);
axis image; axis off;
hold on;
contour(V, label, 'Color', [1 0.2 0.1]);
hold off;

%% 5th Example: restricted domain

% load mask
P = imread('sala.png');
P = rescale(P, 0, 1);
n = size(P,1);
m = size(P,2);

% pick points
S = zeros(n,m);
S = int32(S);
ctrl = figure;
imageplot(P);
title('Pick 5 points in the domain');
axis image; axis off;
hold on;
sp = pick_points(5);
hold off;
sp = int32(sp);

label = 1;
for i=1:size(sp,2)
   S(sp(1,i),sp(2,i)) = label;
   label = label+1;
end

% apply mask: IMPORTANT FOR COMPUTATIONS
S(P==0) = -1;       % negative values in seed matrix are considered outside the domain

% tensor field
T = zeros(n,m,2,2);
T(:,:,1,1) = 1;
T(:,:,2,2) = 1;

% compute distance and Voronoi maps
tic
[U,V,dUx,dUy] = anisoVoronoi2Diterative(T,S);   % 5th USE (with mask implicitely in S)
toc

% display Voronoi diagram
clf(ctrl);
imagesc(V);
axis image; axis off;
hold on;
plot_pts(sp, 'r.',10);
hold off;

