% test for the computation of orthocenters

p = 15;
vertex = rand(2,p);
face = compute_delaunay(vertex);

% compute centers
[q,r] = compute_orthocenter(vertex,face);

% display the triangulation and inner centers
clf;
hold on;
plot_graph(triangulation2adjacency(face),vertex);
plot_circle(q,r);
eta = .1; axis([-eta 1+eta -eta 1+eta]);
