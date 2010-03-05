% test for the computation of in circle property

path(path, 'toolbox/');

p = 400;
rand('state', 12345);
vertex = rand(2,p);

% generate a valid triangulation that is not delaunay 
aniso = 1000;
vertex1 = vertex;
vertex1(1,:) = 1*vertex(1,:) - aniso*vertex(2,:);
vertex1(2,:) = 0*vertex(1,:) - 1*vertex(2,:);
face0 = compute_delaunay(vertex1);

options.display_flips = 0;
tic;
[face, flips, flipsinv] = perform_delaunay_flipping(vertex,face0,options);
toc;

clf;
subplot(1,2,1);
plot_graph(triangulation2adjacency(face0),vertex);
subplot(1,2,2);
plot_graph(triangulation2adjacency(face),vertex);

edge = compute_edges(face);

% check if this is the correct triangulation
if not( triangulation_isequal(face, compute_delaunay(vertex)) )
    warning('This is not the Delaunay triangulation.');
end
% check if the succession of flips, applied to face0, return face
face1 = perform_triangle_flipping(face0, flips);
if not( triangulation_isequal(face1,face) )
    warning('Forward flips did not worked');
end
% check if the reversed succession of flips, applied to face, return face0
face2 = perform_triangle_flipping(face, flipsinv);
if not( triangulation_isequal(face2,face0) )
    warning('Backward flip did not worked.');
end


