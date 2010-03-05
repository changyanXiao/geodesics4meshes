%%
% Test for isotropic/anisotropic remeshing of surfaces

add_base_paths;

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

name = 'bunny';
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
epsilon = .1;

if 0
    % plot tensor field
    clf; hold on;
    plot_mesh(vertex,faces, options);
    v = vertex(:,landmarks);
    rho = .02;
    w = v + rho*Umax1(:,landmarks);
    plot3( [v(1,:);w(1,:)], [v(2,:);w(2,:)], [v(3,:);w(3,:)], 'r');
end

switch metric_type
    case 'isotropic'
        A = (Cmax+Cmin)/2;
        Cmax = A; Cmin = A;
    case 'constant'
        Cmax = ones(n,1);
        Cmin = ones(n,1);
    case 'anisotropic'
        % keep the same.
end

T = perform_surftensor_remapping(Umin,Umax,Cmin,Cmax);


%%
% Test of propagation.

% compute seed matrix (with labels)

landmarks = 1;
m = 400; % number of points to seeds
edges = compute_edges(faces);

% initialize the map
[U, V, dUx, dUy, dUz] = perform_Aniso_Eikonal_Solver_mesh(vertex, faces, T, landmarks);


displist = [100 200 300 400 m]; k = 1;
for i=2:m
    progressbar(i,m);
    % farthest points
    [tmp,landmarks(end+1)] = max(U(:));
    % update
    [U, V, dUx, dUy, dUz] = perform_Aniso_Eikonal_Solver_mesh(vertex, faces, T, landmarks);
    % display
    if k<=length(displist) && displist(k)==i   
        V = round(V);       
        % compute the voronoi triangulation
        faces1 = compute_voronoi_triangulation_mesh(V, faces);
        vertex1 = vertex(:,landmarks);
        options.method = 'slow';
        faces1 = perform_faces_reorientation(vertex1,faces1, options);        
        clf;
        options.start_points = landmarks;
        plot_fast_marching_mesh(vertex,faces, perform_hist_eq(U, 'linear'), [], options);
        saveas(gcf, [rep name '-' metric_type '-dist-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep name '-' metric_type '-dist-' num2str(i) '.png'], 'png');        
        clf;
        plot_fast_marching_mesh(vertex,faces, V, [], options);
        saveas(gcf, [rep name '-' metric_type '-vor-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep name '-' metric_type '-vor-' num2str(i) '.png'], 'png');
        clf;
        plot_mesh(vertex1, faces1, options);
        lighting flat; shading faceted;
      	saveas(gcf, [rep name '-' metric_type '-mesh-' num2str(i) '.eps'], 'epsc');
        saveas(gcf, [rep name '-' metric_type '-mesh-' num2str(i) '.png'], 'png');
        %
        k = k+1;
    end
end
