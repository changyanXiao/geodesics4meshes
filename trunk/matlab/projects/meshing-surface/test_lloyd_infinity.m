%% 
% Test for Lloyd re-centering for Linf norm.

rep = ['../../results/lloyd-linf/'];
if not(exist(rep))
    mkdir(rep);
end

save_eps = 0;
save_png = 1;

%%
% Load the mesh.

name = 'elephant-50kv';
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

%%
% Compute metric.

metric_type = 'anisotropic';
T = compute_surface_metric(vertex,faces, metric_type, options);

%%
% Randomized seeds points.

m = 100;
landmarks = randperm(n); landmarks = landmarks(1:m);

%%
% Lloyd iterations.

niter = 6;
for i=1:niter
    landmark_old = landmarks;
    [landmarks,Ubound,Uland,Q,voronoi_edges] = perform_lloyd_linfty(vertex,faces, T, landmarks, options);
    options.voronoi_edges = voronoi_edges;
    options.start_points = landmark_old;
    % display distance to seeds
    options.startp_col = 'r.';
    clf; plot_fast_marching_mesh(vertex,faces, perform_hist_eq(Uland, 'linear'), [], options);
    if save_eps
        saveas(gcf, [rep name '-lloyd-dland-' num2str(i) '.eps'], 'epsc');
    end
    if save_png
        saveas(gcf, [rep name '-lloyd-dland-' num2str(i) '.png'], 'png');
    end
    % display distance to boundary
    options.startp_col = 'b.';
    clf; plot_fast_marching_mesh(vertex,faces, perform_hist_eq(Ubound, 'linear'), [], options);
    if save_eps
        saveas(gcf, [rep name '-lloyd-dbound-' num2str(i) '.eps'], 'epsc');
    end
    if save_png
        saveas(gcf, [rep name '-lloyd-dbound-' num2str(i) '.png'], 'png');
    end
end


