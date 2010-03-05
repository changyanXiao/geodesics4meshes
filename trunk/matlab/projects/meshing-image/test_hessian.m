% test for Hessian based metrics

path(path,'mex/');
path(path,'toolbox/');

if not(exist('mode'))
mode = 'isotropic';
mode = 'anisotropic';
end

rep = ['results/hessian/' mode '/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Load a smooth image.

n = 400;
sigma = 15 * n/200;

options.bound = 'per';
randn('state', 123);
M = randn(n);
if n==400
    M = circshift(M,[-20 110]);
end
M = perform_blurring(M, 120*n/200, options);
M = perform_blurring( double(M>mean(M(:))), sigma);
U = perform_blurring(randn(n), 120*n/200, options);
M = M + .1*U/std(U(:));
M = rescale(M);

%%
% Compute the metric.

% metric L^p
p = Inf;
q = 1/(1/p+1);

% L2 is in [0,100], most in [0,20]
% L1>L2 is in [0,250] (most in [0,50])
epsilon = 1;

H = compute_hessian(M, options) * n^2;
[e1,e2,l1,l2] = perform_tensor_decomp(H, 'abs');
l1 = abs(l1)+epsilon;
l2 = abs(l2)+epsilon;
switch mode
    case 'isotropic'
        T = perform_tensor_decomp(e1,e2,l1+l2,l1+l2);
    case 'anisotropic'
        D = (l1.*l2).^(q-1);
        T = perform_tensor_decomp(e1,e2,l1.*D,l2.*D);
end

%% 
% Test of propagation.

% compute seed matrix (with labels)

vertex = [n;n]/2;
vertex(:,end+1) = round([n/2;n/3]/2);

S = zeros(n);
S(vertex(1,:)+(vertex(2,:)-1)*n) = 1:size(vertex,2); S = int32(S);
[U,V,dUx, dUy] = anisoVoronoi2Diterative(T,S);   % FIRST USE

%% FP sampling
vertex = [[1;1] [1;n] [n;1] [n;n]]; 
m0 = size(vertex,2);
m = 800;
U = zeros(n);
V = int32(zeros(n));
S = int32(zeros(n));
S(vertex(1,:) + (vertex(2,:)-1)*n) = 1:size(vertex,2);

% initialize the map
[U,V,dUx, dUy] = anisoVoronoi2Diterative(T,S);

displist = [50 100 200 400 m]; k = 1;
for i=m0+1:m
    progressbar(i,m);
    % farthest points
    [tmp,I] = max(U(:));
    [x,y] = ind2sub([n n], I);
    vertex(:,end+1) = [x;y];    
    % update
    p = vertex(:,end);
    U(p(1),p(2)) = 0.0; % update distance
    V(p(1),p(2)) = i;   % update voronoi
    S(p(1),p(2)) = i;   % update seeds
    %
    [dUx, dUy] = anisoVoronoi2Diterative(T,S,U,V);   % FIRST USE
    % display
    if k<=length(displist) && displist(k)==i
        options.col = 'b-'; options.ms = 0; options.ps = 1;
        % distance
        clf; hold on;
        imageplot(perform_hist_eq(U,'linear'));
        h = plot(vertex(2,:), vertex(1,:), 'r.'); set(h, 'MarkerSize', 20);
        colormap jet(256);
        saveas(gcf, [rep mode '-dist-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep mode '-dist-' num2str(i) '.png'], 'png');
        % samples
        clf; hold on;
        imageplot(M);
        h = plot(vertex(2,:), vertex(1,:), 'r.'); set(h, 'MarkerSize', 20);
        saveas(gcf, [rep mode '-img-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep mode '-img-' num2str(i) '.png'], 'png');
        % triangulation
        warning off;
        sel = 1:n-1;
        faces = compute_voronoi_triangulation(double(V), vertex);        
        warning on;
        clf; hold on;
        imageplot(M');        
        plot_graph(triangulation2adjacency(faces),vertex, options);
        h = plot([1 n n 1 1], [1 1 n n 1], 'b.-');
        set(h, 'LineWidth', 2);
        saveas(gcf, [rep mode '-tri-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep mode '-tri-' num2str(i) '.png'], 'png');
        % linear interpolation
        options.verb = 0; options.remove_nan = 0;
        v = compute_orthoproj_triangulation(vertex, faces, M, options);
        M1 = griddata_arbitrary(faces,vertex,v,n, options);
        M1(isnan(M1)) = M(isnan(M1));
        options.ms = 0; options.ps = 1;
        clf; hold on;
        imageplot(M1');
        plot_graph(triangulation2adjacency(faces),vertex, options);
        saveas(gcf, [rep mode '-appr-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep mode '-appr-' num2str(i) '.png'], 'png');
        warning off;
        imwrite(rescale(M1), [rep mode '-appr0-' num2str(i) '.png'], 'png');
        warning on;
        %
        k = k+1;
    end
end
