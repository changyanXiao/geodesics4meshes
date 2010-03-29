% test for image approximation

add_base_paths;

default rep = ['../../results/image-approx/'];
if not(exist(rep))
    mkdir(rep);
end

default metric = ['structure'];
default use_lloyd = 0;
default lloyd_refresh = 10;
default m = 800;    % number of vertices
default n = 256;  % size of the image
default epsilon = 1; 
default alpha = 1;
default sigma_structure = 6*n/512;
default displist = [100 200 400 m]; 

%%
% Load an image.

default name = 'boat';
M = load_image(name,n);
M = rescale(sum(M,3));

str = [name '-' metric '-alpha' num2str(alpha)];

%%
% Compute the metric.
options.bound = 'sym';
switch metric
    case 'hessian'
        H = compute_hessian(M, options);
        [e1,e2,l1,l2] = perform_tensor_decomp(H, 'abs');
    case 'structure'
        H = compute_structure_tensor(M,1,sigma_structure);
        [e2,e1,l1,l2] = perform_tensor_decomp(H, 'abs');
end
l1 = (abs(l1)+epsilon).^alpha;
l2 = (abs(l2)+epsilon).^alpha;
T = perform_tensor_decomp(e1,e2,l1,l2);

%% 
% FP sampling

vertex = [[1;1] [1;n] [n;1] [n;n]]; 
m0 = size(vertex,2);
U = zeros(n);
V = int32(zeros(n)); S = int32(zeros(n));
S(vertex(1,:) + (vertex(2,:)-1)*n) = 1:size(vertex,2);
% initialize the map
[U,V,dUx, dUy] = anisoVoronoi2Diterative(T,S);

k = 1;
progressbar(1,m);
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
    [dUx, dUy] = anisoVoronoi2Diterative(T,S,U,V);   % FIRST USE
    if use_lloyd && mod(i, lloyd_refresh)==0
        %% LLOYD ITERATIONS %%
        [Y,X] = meshgrid(1:n,1:n);
        for j=1:size(vertex,2)
            I = find(V==j);
            landmark(:,j) = [mean(X(I));mean(Y(I))];
        end
        % recompute distance map
        U = zeros(n);   V = int32(zeros(n)); S = int32(zeros(n));
        S(vertex(1,:) + (vertex(2,:)-1)*n) = 1:size(vertex,2);
        [U,V,dUx, dUy] = anisoVoronoi2Diterative(T,S);
    end
    % display
    if k<=length(displist) && displist(k)==i
        options.col = 'b-'; options.ms = 0; options.ps = 1;
        % distance
        if 0
        clf; hold on;
        imageplot(perform_hist_eq(U,'linear'));
        h = plot(vertex(2,:), vertex(1,:), 'r.'); set(h, 'MarkerSize', 20);
        colormap jet(256);
        saveas(gcf, [rep str '-dist-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep str '-dist-' num2str(i) '.png'], 'png');
        % samples
        clf; hold on;
        imageplot(M);
        h = plot(vertex(2,:), vertex(1,:), 'r.'); set(h, 'MarkerSize', 20);
        saveas(gcf, [rep str '-img-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep str '-img-' num2str(i) '.png'], 'png');
        end
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
        saveas(gcf, [rep str '-tri-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep str '-tri-' num2str(i) '.png'], 'png');
        % linear interpolation
        options.verb = 0; options.remove_nan = 0;
        v = compute_orthoproj_triangulation(vertex, faces, M, options);
        M1 = griddata_arbitrary(faces,vertex,v,n, options);
        M1(isnan(M1)) = M(isnan(M1));
        if 0
        options.ms = 0; options.ps = 1;
        clf; hold on;
        imageplot(M1');
        plot_graph(triangulation2adjacency(faces),vertex, options);
        saveas(gcf, [rep mode '-appr-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep mode '-appr-' num2str(i) '.png'], 'png');
        warning off;
        imwrite(rescale(M1), [rep mode '-appr0-' num2str(i) '.png'], 'png');
        warning on;
        end
        %
        k = k+1;
    end
end

%%
% compute the approximation on the triangulation.
warning off;
faces = compute_voronoi_triangulation(double(V), vertex);
warning on;
options.verb = 0; options.remove_nan = 0;
v = compute_orthoproj_triangulation(vertex, faces, M, options);
M1 = griddata_arbitrary(faces,vertex,v,n, options);
M1(isnan(M1)) = M(isnan(M1));
