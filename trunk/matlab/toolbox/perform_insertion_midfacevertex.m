function [faces, vertex, fring] = perform_insertion_midfacevertex(faces, vertex, fring, i)

% perform_insertion_midfacevertex - add a vertex in the middle of a face
%
%   [faces, vertex, fring] = perform_insertion_midfacevertex(faces, vertex, fring, i)
%
%   add vertex(:,end+1) in the middle of faces(:,i), and update the
%   faces/vertex/fring matrices.
%
%   Copyright (c) 2010 Gabriel Peyre

n = size(vertex,2);
m = size(faces,2);

f = faces(:,i);
%% add the vertex
vertex(:,n+1) = mean(vertex(:,f), 2);
%% add the three faces
faces(:,i) = [f(1);n+1;f(3)];
faces(:,m+1) = [f(2);n+1;f(1)];
faces(:,m+2) = [f(3);n+1;f(2)];
%% update fring
r = fring(:,i);
% the new
fring(:,i) = [m+2; r(2); m+1];
fring(:,m+1) = [i; r(3); m+2];
fring(:,m+2) = [m+1; r(1); i];
% the olds
if r(1)>0
    I = find(fring(:,r(1))==i);
    fring(I,r(1)) = m+2;
end
if r(3)>0
    I = find(fring(:,r(3))==i);
    fring(I,r(3)) = m+1;
end


