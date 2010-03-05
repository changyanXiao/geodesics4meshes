function face_voronoi = compute_voronoi_triangulation_mesh(Q, faces)

% compute_voronoi_triangulation_mesh - compute a triangulation
%
%   face_voronoi = compute_voronoi_triangulation_mesh(Q, vertex, faces);
%
%   Q is a Voronoi partition function, computed using
%   perform_fast_marching_mesh.
%
%   Copyright (c) 2006 Gabriel Peyre

if size(faces,1)>size(faces,2)
    faces = faces';
end

V = Q(faces);
V = sort(V,1);
V = unique(V', 'rows')';
% V = V( prod(V,2)>0 ,:);

% remove inf
V = V(:,sum(V,1)~=Inf);

d = (V(1,:)~=V(2,:)) + (V(2,:)~=V(3,:));

I = find(d==2); I = sort(I);

face_voronoi = V(:,I);

return;

%% OLD %%

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end


w = V(:,I); w = sort(w(:)); % index that are in the triangulation

nverts = size(vertex,2);
z = zeros(nverts,1);
z(w) = (1:length(w))';

faces_voronoi = z(V( :,I ));
% vertex_voronoi = vertex(:,w);
