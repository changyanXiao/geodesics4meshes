%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_aniso_contour_completion.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path, 'data/');
name = 'curves_b';
M = load_image(name);
M = rescale(M);
n = size(M,1);
m = size(M,2);
epsilon = 0;
%% compute structure tensor around the curve
alpha = 1;
sigma1 = 1;   % for gradient computation
sigma2 = 6;     % for filtering the tensor field
%options.use_anisotropic = 0;
T = compute_structure_tensor(M, sigma1, sigma2);
%T = compute_hessian_gaussian(M, 5, 6);
%T = compute_hessian(M);
%[e1,e2,l1,l2] = perform_tensor_decomp(T);
%T = perform_tensor_recomp(e2, e1, ((l1+epsilon).^alpha), ((l2+epsilon).^alpha));
%T = perform_tensor_recomp(e1, e2, (abs(l1)+epsilon).^alpha,(abs(l2)+epsilon).^alpha);
%T = perform_tensor_recomp(e1,e2,l1*0+1,l2*0);

%plot_tensor_field_color(T);
%% mask
I = zeros(n);
%I(M<0.9) = 1;
%[e1,e2,l1,l2] = perform_tensor_decomp(T);
%l1(I==0)=0;
%l2(I==0)=0;
%l1(I==1)=1;
%l2(I==1)=0;
%T1 = perform_tensor_recomp(e1,e2,l1,l2);

%plot_tensor_field_color(T1);
%% tensor densification
T1 = tensor22_to_tensor3(T);
options.niter = 600;
options.sigma = [16,1.3];
T1 = perform_tensor_densification(T1, I, options);
T1 = tensor3_to_tensor22(T1);
figure;
[e1,e2,l1,l2] = perform_tensor_decomp(T1);
T = perform_tensor_recomp(e1,e2,l1,l2);

%plot_tensor_field(T);
%% get points
P = pick_points(M, 6);
%save('points','P');
%load 'points'
%% Reconstruction
% SP=saddle points, E=edges
[U,SP,V,E] = AnisoContourCompletion2D(T1, P);
V=V+1;
E=E+1;
Uh = perform_histogram_equalization(U, linspace(0,255,size(M,1)*size(M,2)));
W = double(V);
%% Distance map and combinatorial edges
%faces=compute_voronoi_triangulation(W,P);
% plot Voronoi label map
figure;
imagesc(Uh);
hold on;
axis image; axis off;
plot_edges(E, P, 'w', 4);
plot_pts(P,'w.',24);
colormap jet(255);
hold off;
%% Voronoi diagram
I=V;
k = mmax(V)*1.3;
for i=2:n-1
    for j=2:n-1
        if V(i,j)~=V(i+1,j)
            I(i,j)=k;
            I(i+1,j)=k;
        end
        if V(i,j)~=V(i-1,j)
            I(i,j)=k;
            I(i-1,j)=k;
        end
        if V(i,j)~=V(i-1,j-1)
            I(i,j)=k;
            I(i-1,j-1)=k;
        end
        if V(i,j)~=V(i-1,j+1)
            I(i,j)=k;
            I(i-1,j+1)=k;
        end
        if V(i,j)~=V(i+1,j-1)
            I(i,j)=k;
            I(i+1,j-1)=k;
        end
        if V(i,j)~=V(i,j-1)
            I(i,j)=k;
            I(i,j-1)=k;
        end
        if V(i,j)~=V(i,j+1)
            I(i,j)=k;
            I(i,j+1)=k;
        end
        if V(i,j)~=V(i+1,j+1)
            I(i,j)=k;
            I(i+1,j+1)=k;
        end
    end
end
figure;
imagesc(I);
axis image; axis off;
colormap hot(256);
hold on;
plot_edges(edges, P, 'k', 4);
plot_pts(P,'w.',26);
hold off;
%% Geodesics
% compute and plot paths
npaths=size(SP,2);
%options.trim_path=0;
paths={};
figure;
imagesc(M);
axis image; axis off;
colormap jet(256);
hold on;
for i=1:npaths
    %A=U;
    %B=find(V~=V(SP(1,i),SP(2,i)));
    %A(B)=Inf;
    paths{i} = compute_geodesic(U,SP(:,i));
    h = plot( paths{i}(2,:), paths{i}(1,:), 'k' );
    set(h, 'LineWidth', 3);
end
% plot saddle_points
j=1;
stop = npaths/2;
for i=1:stop
    h = plot((SP(2,j)+SP(2,j+1))/2.,(SP(1,j)+SP(1,j+1))/2., 'k+');
    set(h, 'MarkerSize', 12);
    set(h, 'LineWidth', 4);
    j=j+2;
end
plot_pts(P,'w.',24);
%saveas(gcf, [rep name '-recons.eps'], 'psc2');
hold off;
