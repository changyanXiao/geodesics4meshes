%%
% Test for isotropic/anisotropic remeshing of surfaces

add_base_paths;
path(path, '../../data/meshes-large/');


save_eps = 0;
save_png = 1;

%if not(exist('metric_type'))
metric_type = 'isotropic';
metric_type = 'constant';
metric_type = 'anisotropic';
%end

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

T = compute_surface_metric(vertex,faces, metric_type, options);

if 0
    u = Umin(:,1);
    v = Umax(:,1);
    T = perform_surftensor_remapping([u v],[v u],[.1 1],[.1 1]);
    
    T1 = reshape(T, [6 1 size(T,2)]);
    H = zeros(3,3,size(T,2));
    H(1,1,:) = T1(1,1,:);
    H(2,2,:) = T1(2,1,:);
    H(3,3,:) = T1(3,1,:);
    H(1,2,:) = T1(4,1,:);
    H(2,3,:) = T1(5,1,:);
    H(1,3,:) = T1(6,1,:);
    H(2,1,:) = H(1,2,:);
    H(3,1,:) = H(1,3,:);
    H(3,2,:) = H(2,3,:);
end

if 0
    clf;
    options.face_vertex_color = rescale(Cmin(:));
    plot_mesh(vertex,faces, options);
    colormap jet(256);
end

%%
% Test of propagation.

% compute seed matrix (with labels)

landmarks = 1;

% initialize the map
Calls = 0;
options.doUpdate = true(n,1);
[U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks);
ULoc(:,1) = U;

% smoothing matrix, useful to extend the voronoi
W = triangulation2adjacency(faces) + speye(n);
% size of the extension for the voronois
Vext = 2;

% rate for the Lloyd iterations
recentering_rate = 20;

m = 6400; % number of points to seeds
displist = 100 * 2.^(1:log2(m/100));
k = 1;
for i=2:m
    progressbar(i,m);
    [tmp,landmarks(end+1)] = max(U(:));
    
    if 0
        % compute triple point
        [flist,lambda] = compute_double_points(vertex,face, landmarks, Uloc);
        if isempty(flist)
            % take point with farthest distance
            d = sum(U(flist).*lambda, 3);
            [tmp,i0] = max(d); f = flist(:,i0); l = lambda(:,i0);
            % add the point to the mesh
            vertex(:,end+1) = vertex(:,f).*l;
            n = n+1;
        else
            % farthest points
            [tmp,landmarks(end+1)] = max(U(:));
        end
    end
    
    % update of the distance map
    options.doUpdate = true(n,1);
    options.U_ini = U;
    options.V_ini = V;
    [U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks(end), options);

    if 0
        % compute a slightly enlarged voronoi region
        v = double(V==i);
        for iext=1:Vext
            v = W*v;
        end
        v = v>0;
        % perform propagation on the enlarged region
        options.doUpdate = v;
        options.U_ini = []; % zeros(n,1);
        options.V_ini = []; % ones(n,1);
        [ULoc(:,i), Vloc, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks(end), options);
        ULoc(v==0,i) = Inf;
    end
    
    % Perform re-centering
    if mod(i,recentering_rate)==0
        [landmarks,U,Uland,V,voronoi_edges] = perform_lloyd_linfty(vertex,faces, T, landmarks, options);
        V = int32(V);
    end
    
    % display
    if k<=length(displist) && displist(k)==i
        Va = round(double(V));
        % compute the voronoi triangulation
        faces1 = compute_voronoi_triangulation_mesh(Va, faces);
        vertex1 = vertex(:,landmarks);
        options.method = 'slow';
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
            % display voronoi with there boundaries.
            options.Dlist = ULoc;
            [Q, DQ, voronoi_edges, edges_id, lambda] = compute_voronoi_mesh(vertex,faces, landmarks, options);
            clf;
            options.voronoi_edges = voronoi_edges;
            plot_fast_marching_mesh(vertex,faces, Q(:,1), [], options);
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
        % lighting flat;
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
    
    if i==400 && strcmp(metric_type,'constant')
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
