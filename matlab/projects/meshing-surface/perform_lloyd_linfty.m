function  [landmarks1,Ubound,Uland,Q,voronoi_edges, Calls] = perform_lloyd_linfty(Calls, vertex,faces, T, landmarks, options)

% perform_lloyd_linfty - perform lloyd recentering
%
%   [landmarks1,U,Q,voronoi_edges] = perform_lloyd_linfty(vertex,faces, T, landmarks, options)
%
%   Replaces each landmarks(i) by the index landmarks1(i) which is the
%   points in voronoi(i) that is the farthest away from the boundary of
%   voronoi(i).
%
%   Useful for display purpose: 
%       U is the distance to the boundary
%       Q is the voronoi segmentation
%       voronoi_edges is the set of edges of the boundary of the segmentation.
%
%   Set options.lloyd_niter to make more than 1 iteration.
%
%   Copyright (c) Gabriel Peyre 2010.

options.null = 0;
lloyd_niter = getoptions(options, 'lloyd_niter', 1);
Vext = getoptions(options, 'vornoi_extension', 3);

if lloyd_niter>1
    options.lloyd_niter = 1;
    for i=1:lloyd_niter
        [landmarks,U,Q,voronoi_edges] = perform_lloyd_linfty(vertex,faces, T, landmarks, options);
    end
    landmarks1 = landmarks;
    return;
end

n = size(vertex,2);
m = length(landmarks);

% Smoothing matrix, useful to extend the voronoi.

W = triangulation2adjacency(faces) + speye(n);

% Compute pixelic Voronois.

options.doUpdate = true(n,1);
[U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarks);
 
% Compute distance to seeds with overlap
ULoc = zeros(n,m);
for i=1:m
    [ULoc(:,i), Calls] = perform_geodesic_computation_extended(Calls, vertex, faces, W, T, U, V, i, Vext);
end

Uland = min(ULoc,[],2);

% Compute Voronoi boundaries.

options.Dlist = ULoc;
options.Dland = Uland;

[Q, DQ, voronoi_edges, edges_id, lambda] = compute_voronoi_aniso_mesh(vertex,faces, landmarks, options);
Q = Q(:,1);

% Perform initialization along the boundaries.

% compute the length of each edge
ne = size(edges_id,2);
le = zeros(ne,1);
% mean tensor along the edges
t = ( T(:,edges_id(1,:)) + T(:,edges_id(2,:)) )/2;
% voronoi edges in 3D
v = vertex(:,edges_id(1,:)) - vertex(:,edges_id(2,:));
% length according to metric
le = compute_tensor_distance(v, t);
% initial distance by 
nbr = zeros(n,1);
d = zeros(n,1);
for i=1:ne
    e = edges_id(:,i);
    nbr(e) = nbr(e)+1;
    d(e(1)) = d(e(1)) + lambda(i) * le(i);
    d(e(2)) = d(e(2)) + (1-lambda(i)) * le(i);
end
% seeds for the propagation
landmarksV = find(nbr>0);
distV = d(nbr>0) ./ nbr(nbr>0);


% Perform FM from the boundaries.

q = length(landmarksV);
options.doUpdate = true(n,1);
options.U_ini = []; % zeros(n,1);
options.V_ini = []; % ones(n,1);
options.U_ini_seeds = distV(1:q);
options.V_ini_seeds = int16(distV(1:q));%just for allocations: this value is useless take whatever but at the right size
[Ubound, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, landmarksV(1:q), options);
options.U_ini_seeds = [];

% Re-center by farthest points inside each voronoi region.

landmarks1 = landmarks;
for i=1:m
    I = find(Q==i);
    [tmp,i0] = max(Ubound(I));
    landmarks1(i) = I(i0);
end

Calls = Calls+1;
return;