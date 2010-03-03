%%
path(path, 'data/');
%name = 'mountain';
name = 'eagle';
%n=600;
M = load_image(name);
M = rescale(M);
n = size(M,1);
m = size(M,2);
%m=n;
alpha = 1.0;
epsilon = 0.000001;
%%
% compute structure tensor around the curve
sigma1 = 1.5;   % for gradient computation
sigma2 = 2;     % for filtering the tensor field
options.use_anisotropic = 0;
T = compute_structure_tensor(M, sigma1, sigma2);
[e1,e2,l1,l2] = perform_tensor_decomp(T);
T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);
%%
nb_pts_to_insert = 50;

[U,dUx,dUy,V,Pts,faces] = FPSCAniso2D([1;1], T, nb_pts_to_insert, 0);

Uf = perform_histogram_equalization(U, linspace(0,1,size(M,1)*size(M,2)));
V=V+1;
faces=faces+1;
%%
figure;
imagesc(V);
hold on;
axis image; axis off;
colormap jet(256);
edges = compute_edges(faces);
plot_edges(edges, Pts, 'k', 3);
plot_pts(Pts,'k.',20);
hold off;
