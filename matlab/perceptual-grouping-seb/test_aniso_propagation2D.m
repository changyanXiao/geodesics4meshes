%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_aniso_propagation2D.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path, 'data/');
name = 'varying-2d-texture';
M = load_image(name);
M = rescale(M);
alpha = 1;
n = size(M,1);
m = size(M,2);
epsilon = 0.000001;
%%
% compute structure tensor around the curve
sigma1 = 1.5;   % for gradient computation
sigma2 = 4;     % for filtering the tensor field
T = compute_structure_tensor(M, sigma1, sigma2);
[e1,e2,l1,l2] = perform_tensor_decomp(T);
T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);
%%
% identity
%T = ones(n,n,2,2);
%T(:,:,1,1) = 1;
%T(:,:,2,2) = 1;
%T(:,:,1,2) = 0;
%T(:,:,2,1) = 0;
%%
% get seeds
nb_seeds = 10;
P = pick_points(M, nb_seeds);
%%
% propagations
eps_update = 0;  % positif et proche de zero (lie a la convergence de la version resursive),
                 % plus il est grand et plus la version recursive est
                 % proche de la non recursive.
[U1,V1] = RecursiveAnisoPropagation2D(T, P, 0);
[U2,V2] = AnisoPropagation2D(T, P);
V=V+1;
%%
% plot results
figure;
U1 = perform_histogram_equalization(U1, linspace(0,1,size(M,1)*size(M,2)));
V1 = double(V1);
imageplot({U1, V1},{'recursive propagation', 'Voronoi'});
hold on;
axis image; axis off;
colormap jet(256);
plot_pts(P,'k.',28);
hold off;
figure;
U2 = perform_histogram_equalization(U2, linspace(0,255,size(M,1)*size(M,2)));
V2 = double(V2);
imageplot({U2, V2},{'propagation', 'Voronoi'});
hold on;
axis image; axis off;
colormap jet(256);
plot_pts(P,'k.',28);
hold off;
