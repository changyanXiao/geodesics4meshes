%%
% test for approximate Lloyd algorithm, using Euclidean distance for the
% re-centering.

rep = ['../../results/lloyd-approx/'];
if not(exist(rep))
    mkdir(rep);
end

save_eps = 0;
save_png = 1;
add_base_paths;

%%
% Load the mesh.

name = 'elephant-50kv';
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

%%
% Compute metric.

metric_type = 'anisotropic';
metric_type = 'constant';
T = compute_surface_metric(vertex,faces, metric_type, options);

%%
% Randomized seeds points.


m = 300;
landmarks = randperm(n); landmarks = landmarks(1:m);


I = find(vertex(2,:)>median(vertex(2,:)));
landmarks = randperm(length(I));
landmarks = I(landmarks(1:m));

%%
% Lloyd iterations.
Calls = 0;
niter = 128;
displist = [1 2 4 8 16 32 64 niter]; k =1;
for i=1:niter
    
    %%
    % Compute Voronoi partition
    
    options.doUpdate = true(n,1);
    [U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks);
    
    if i==displist(k)
        %        [Q, DQ, voronoi_edges, edges_id, lambda] = compute_voronoi_mesh(vertex,faces, landmarks);
        %        options.voronoi_edges = voronoi_edges;
        clf;
        options.startp_col = 'k.';
        options.v_edge_color = 'm-';
        options.start_points = landmarks;
        plot_fast_marching_mesh(vertex,faces, double(V), [], options);
        k = k+1;
        if save_eps
            saveas(gcf, [rep name '-lloyd-' num2str(i) '.eps'], 'epsc');
        end
        if save_png
            print2im([rep name '-lloyd-' num2str(i) '.png'], '-alpha');
        end
    end
    
    %%
    % perform re-centering
    
    for j=1:m
        I = find(V==j);
        W = vertex(:,I);
        % centeroid
        w = mean(W,2);
        % distance to centroid
        d = sum( (W - repmat(w, [1 length(I)])).^2 );
        % nearest point
        [tmp,nn] = min(d);
        landmarks(j) = I(nn);
    end
end