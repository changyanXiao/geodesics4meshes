% test for in-circle test

n = 100;
vertex = randn(2,n);
face = compute_delaunay(vertex);

% every edge should pass the test
T = check_incircle_edge(vertex,face);

% should be 0
sum(T~=1)