%%
% Test for isotropic/anisotropic remeshing of surfaces

add_base_paths;
path(path, '../../data/meshes-large/');

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

options.verb = 0;
[Umin,Umax,Cmin,Cmax] = compute_curvature(vertex,faces,options);
Cmin = abs(Cmin); Cmax = abs(Cmax);
% ensure Cmin<Cmax
I = find(Cmin>Cmax);
Cmax1=Cmax; Cmax(I) = Cmin(I); Cmin(I) = Cmax1(I);
Umax1=Umax; Umax(:,I) = Umin(:,I); Umin(:,I) = Umax1(:,I);
% ensure strong convexity
%epsilon = .000;
%Cmin = Cmin+epsilon; Cmax = Cmax + epsilon;

A = Cmax./Cmin;
E = Cmax+Cmin;

switch metric_type
    case 'isotropic'
        A = A*0+1;
        E = perform_hist_eq(A, 'linear');
        E = rescale(E,1,10);
    case 'constant'
        A = A*0+1;
        E = E*0+1; 
    case 'anisotropic'
%        E = E*0+1;
%        A = A*0+100;
end

Cmin = E./(1+A);
Cmax = Cmin.*A;

T = perform_surftensor_remapping(Umin,Umax,Cmin,Cmax);

% Cmax = rescale(-vertex(1,:)>0, .01,1); Cmin = Cmax;
% T = perform_surftensor_remapping(Umin,Umax,Cmin,Cmax);

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
        
m = 6400; % number of points to seeds
displist = 100 * 2.^(1:log2(m/100));
k = 1;
for i=2:m
    progressbar(i,m);
    % farthest points
    [tmp,landmarks(end+1)] = max(U(:));
    % update everywhere
    options.doUpdate = true(n,1);
    options.U_ini = U;
    options.V_ini = V;
    [U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks(end), options);
    % compute a slightly enlarge voronoi region
    v = double(V==i);
    for iext=1:Vext
        v = W*v;
    end
    v = v>0;
    options.doUpdate = v;
    options.U_ini = []; % zeros(n,1);
    options.V_ini = []; % ones(n,1);
    [ULoc(:,i), Vloc, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks(end), options);
    ULoc(v==0,i) = Inf;    
    % compute Voronoi diagram in sub-pixelic way
    options.Dlist = ULoc;
    [Q, DQ, voronoi_edges, edges_id, lambda] = compute_voronoi_mesh(vertex,faces, landmarks, options);
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
            saveas(gcf, [rep name '-' metric_type '-dist-' num2str(i) '.eps'], 'epsc');
            saveas(gcf, [rep name '-' metric_type '-dist-' num2str(i) '.png'], 'png');
            % display voronoi with there boundaries.
            clf;
            options.voronoi_edges = voronoi_edges;
            plot_fast_marching_mesh(vertex,faces, Q(:,1), [], options);
            saveas(gcf, [rep name '-' metric_type '-vor-' num2str(i) '.eps'], 'epsc');
            saveas(gcf, [rep name '-' metric_type '-vor-' num2str(i) '.png'], 'png');
            options.voronoi_edges = [];
        end
        clf;
        options.face_vertex_color = [];
        plot_mesh(vertex1, faces1, options);
        % lighting flat; 
        shading faceted;
      	saveas(gcf, [rep name '-' metric_type '-mesh-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep name '-' metric_type '-mesh-' num2str(i) '.png'], 'png');
        %
        k = k+1;
    end
    
    if i==400 && strcmp(metric_type,'constant')
        clf;
        plot_surface_tensor(vertex,faces,{Umax Umin},landmarks, options);
        saveas(gcf, [rep name '-tensor.eps'], 'epsc');
        saveas(gcf, [rep name '-tensor.png'], 'png');
    end
end
