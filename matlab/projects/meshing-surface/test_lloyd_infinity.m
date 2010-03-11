%% 
% Test for Lloyd re-centering for Linf norm.
clear all; close all; clc;
rep = ['../results/lloyd-linf/'];
if not(exist(rep))
    mkdir(rep);
end

save_eps = 0;
save_png = 1;
%%
library_directory = '../mex';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
meshes_directory   = '../meshes/';
addpath(meshes_directory);
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
Calls = 0;
niter = 10;
for i=1:niter
    landmark_old = landmarks;
    tic
    [landmarks,Ubound,Uland,Q,voronoi_edges, Calls] = perform_lloyd_linfty(Calls, vertex,faces, T, landmarks, options);
    toc
    options.voronoi_edges = voronoi_edges;
    options.start_points = landmark_old;
    % display distance to seeds
    options.startp_col = 'r.';
    options.v_edge_color = 'm-';
    clf; plot_fast_marching_mesh(vertex,faces, perform_hist_eq(Uland, 'linear'), [], options);
    if save_eps
        saveas(gcf, [rep name '-lloyd-dland-' num2str(i) '.eps'], 'epsc');
    end
    if save_png
        print2im([rep name '-lloyd-dland-' num2str(i) '.png'], '-alpha');
    end
    % display distance to boundary
    options.startp_col = 'b.';
    options.v_edge_color = 'r-';
    clf; plot_fast_marching_mesh(vertex,faces, perform_hist_eq(Ubound, 'linear'), [], options);
    if save_eps
        saveas(gcf, [rep name '-lloyd-dbound-' num2str(i) '.eps'], 'epsc');
    end
    if save_png
        print2im([rep name '-lloyd-dbound-' num2str(i) '.png'], '-alpha');
    end
end