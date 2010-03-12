% test for point insertion inside a surface.


add_base_paths;

name = 'beetle';
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

clf;
plot_mesh(vertex, faces, options);
shading faceted;

fring = compute_face_ring_enhanced(faces);

%% Add a new vertex in the middle of a face.

for i=1:1000
    i = floor(rand*size(vertex,2))+1;
    [faces, vertex, fring] = perform_insertion_midfacevertex(faces, vertex, fring, i);
%    I = find(sum(abs(fring-fring1))>0);
end

clf;
plot_mesh(vertex, faces, options);

%% 
% The final test

fring1 = compute_face_ring_enhanced(faces);
disp(['This should be 0: ' num2str(norm(fring-fring1, 'fro')) '.']);

shading faceted;