function fring = compute_face_ring_enhanced(faces)

% compute_face_ring_enhanced - same as compute_face_ring but better
%
%   fring = compute_face_ring_enhanced(vertex,faces);
%
%   r = fring(:,i) = [f1;f2;f3] is the set of faces ajacent to face i, 
%   with the convention that f2=-1 and/or f3=-1 if the face is on the
%   boundary of the mesh.
%
%   The adjacent faces are order so that if faces(:,i)=[v1;v2;v3], then
%       * [v2;v3] is an edge of f1
%       * [v1;v3] is an edge of f2
%       * [v1;v2] is an edge of f3
%
%   This allows for a simple traversal of the faces.
%
%   Copyright (c) 2010 Gabriel Peyre

m = size(faces,2);
A = compute_edge_face_ring(faces);
[i,j,s1] = find(A);     % direct link
[i,j,s2] = find(A');    % reverse link
I = find(i<j);
s1 = s1(I); s2 = s2(I);
% build face ring
fring = zeros(3, m)-1;
cnt = ones(m,1);
for k=1:length(s1)
    if s1(k)>0 && s2(k)>0
        fring(cnt(s1(k)),s1(k)) = s2(k);
        cnt(s1(k)) = cnt(s1(k))+1;
        fring(cnt(s2(k)),s2(k)) = s1(k);
        cnt(s2(k)) = cnt(s2(k))+1;
    end
end

% re-ordering of the faces
for i=1:m
    f = fring(:,i);
    f1 = -ones(3,1);
    for j=1:3
        if f(j)>0
            % common edge
            e = sort(intersect(faces(:,f(j)),faces(:,i)));
            if norm(e-sort(faces(2:3,i)))==0
                f1(1) = f(j);
            elseif norm(e-sort(faces([1 3],i)))==0
                f1(2) = f(j);
            elseif norm(e-sort(faces(1:2,i)))==0
                f1(3) = f(j);
            else
                error('Problem with the triangulation');
            end
        end
    end
    fring(:,i) = f1;        
end