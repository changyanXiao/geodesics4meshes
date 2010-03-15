%%
% Test for isotropic/anisotropic remeshing of surfaces
% using pixelic precision (not the best).

add_base_paths;
path(path, '../../data/meshes-large/');

save_eps = 0;
save_png = 1;

if not(exist('metric_type'))
metric_type = 'isotropic';
metric_type = 'constant';
metric_type = 'anisotropic';
metric_type = 'isotropic-boost';
metric_type = 'anisotropic-boost';
end

rep = ['../../results/meshing-surface/' metric_type '/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Load a mesh.

name = 'vase-lion';
name = 'bunny';
name = 'kitten';
name = 'rocker-arm';
name = 'screwdriver';
name = 'fertility';
name = 'elephant-50kv';
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

%%
% Re-center
vertex = vertex - repmat(mean(vertex,2), [1 n]);
vertex = vertex/max(abs(vertex(:)));

%%
% Compute the metric.

[T,Umin,Umax] = compute_surface_metric(vertex,faces, metric_type, options);

%%
% Test of propagation.

% compute seed matrix (with labels)
landmarks = 1;

% initialize the map
Calls = 0;
options.doUpdate = true(n,1);
[U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks);

% rate for the Lloyd iterations
recentering_rate = 100;

m = 6400/2; % number of points to seeds
displist = 100 * 2.^(2:log2(m/100));
k = 1;
for i=2:m
    progressbar(i,m);

    % farthest points
    [tmp,landmarks(end+1)] = max(U(:));
    
    % update of the distance map
    options.doUpdate = true(n,1);
    options.U_ini = U;
    options.V_ini = V;
    [U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks(end), options);    
    if mod(i,recentering_rate)==0
        % Perform re-centering
        [landmarks,Ubound,Uland,Vold,voronoi_edges] = perform_lloyd_linfty(Calls, vertex,faces, T, landmarks, options);
        % recompute distance map
        options.doUpdate = true(n,1);
        options.U_ini = [];
        options.V_ini = [];
        [U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks);
    end    
    % display
    if k<=length(displist) && displist(k)==i
        Va = round(double(V));
        % compute the voronoi triangulation
        faces1 = compute_voronoi_triangulation_mesh(Va, faces);
        vertex1 = vertex(:,landmarks);
        options.method = 'slow';
        options.verb = 0;
        faces1 = perform_faces_reorientation(vertex1,faces1, options);
        if i<=800
            clf;
            options.start_points = landmarks;
            plot_fast_marching_mesh(vertex,faces, perform_hist_eq(U, 'linear'), [], options);
            if save_eps
                saveas(gcf, [rep name '-' metric_type '-dist-' num2str(i) '.eps'], 'epsc');
            end
            if save_png
                saveas(gcf, [rep name '-' metric_type '-dist-' num2str(i) '.png'], 'png');
            end
            
            if 0
            clf;
            plot_fast_marching_mesh(vertex,faces, V, [], options);
            if save_eps
                saveas(gcf, [rep name '-' metric_type '-vor-' num2str(i) '.eps'], 'epsc');
            end
            if save_png
                saveas(gcf, [rep name '-' metric_type '-vor-' num2str(i) '.png'], 'png');
            end
            options.voronoi_edges = [];
            end
        end
        clf;
        options.face_vertex_color = [];
        plot_mesh(vertex1, faces1, options);
        shading faceted;
        if save_eps
            saveas(gcf, [rep name '-' metric_type '-mesh-' num2str(i) '.eps'], 'epsc');
        end
        if save_png
            saveas(gcf, [rep name '-' metric_type '-mesh-' num2str(i) '.png'], 'png');
        end
        %
        k = k+1;
    end
    
    if i==800 && strcmp(metric_type,'constant')
        clf;
        plot_surface_tensor(vertex,faces,{Umax Umin},landmarks, options);
        if save_eps
            saveas(gcf, [rep name '-tensor.eps'], 'epsc');
        end
        if save_png
            saveas(gcf, [rep name '-tensor.png'], 'png');
        end
    end
end
