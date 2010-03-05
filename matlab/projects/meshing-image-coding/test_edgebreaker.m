%% test for edgebreaker

%% compute a random delaunay triangulation
n = 300;
vertex = rand(2,n);
face = compute_delaunay(vertex);

%% display the triangulation
options.verb = 0;
bound = compute_boundary(face, options);
clf;
hold on;
plot_mesh(vertex,face);
h = plot(vertex(1,bound), vertex(2,bound), 'r.-');
set(h, 'LineWidth', 2);

%% Options for EB
% mesh with boundary
options.patch = 1;

%% First pass - ensure that vertex and face are in the correct EB order %%
% [stream,vertexC] = perform_edgebreaker(face, vertex, +1, options);
% [face,vertex] = perform_edgebreaker(stream, vertexC, -1, options);

%% CODE %%
% note : vertexC is useless, just 0's
[stream,vertexC] = perform_edgebreaker(face, vertex, +1, options);
%% DECODE %%
[face1,vertex1] = perform_edgebreaker(stream, vertexC, -1, options);

% display
clf;
plot_mesh(vertex1,face1);


if 0
    ov_filename = 'OVTable_Cow.ov';
    ov_filename_out = 'OVTable_Cow_decomp.ov';
    output_dir = '.';
    EBdir = './';
    surf_type = 'MANIFOLD';
    StartCorner = 1;
    format = 'ASKII';

    system([EBdir 'EBCompression ' ov_filename ' ' surf_type ' ' num2str(StartCorner)  ' ' format ]);
    system([EBdir 'EBDecompression ' output_dir ' ' ov_filename_out ' ' format]);
end
