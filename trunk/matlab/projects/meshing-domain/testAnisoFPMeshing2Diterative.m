% clear all; close all; clc;
%--------------------------------------------------------------------------
add_base_paths;

%%
% load mask
P = imread('sala.png');
n = size(P,1);
m = size(P,2);

n = 300;
M = rescale(load_image('sala',n));
M = perform_blurring(M,5)>.5;
M = conv2(double(M), ones(3), 'same')>0;
P = double(M);


n = size(P,1);
m = size(P,2);


% SE = strel('square',3);
% P = imopen(P,SE);    % IMPORTANT : le bord du domaine doit pouvoir s'extraire en 4-conexité

P = conv2(double(P), ones(3), 'same')>0;
% P = conv2(double(P==0), ones(3), 'same')==0;

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
nb_seeds = 300;
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

% convert voronoi to color
[a,b,W] = unique(V(:));
W = reshape(W, size(V));
c = jet(256);
I = floor(255*rescale(W))+1;
Z = c(I(:),:); Z = reshape(Z, [size(W) 3]);
I = find(P==0);
Z([I I+prod(size(P)) I+2*prod(size(P))]) = 1;

imageplot(Z);
colormap jet(256);
hold on;
plot_pts(sp,'k.',20);
hold off;


return;

imagesc(U);
axis image; axis off;
colormap jet(256);
hold on;
plot_pts(sp, 'r.',10);
hold off;

