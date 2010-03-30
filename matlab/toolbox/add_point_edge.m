function [vertex, connectivity, modified] = add_point_edge(vertex, connectivity, edge, coords)

% add_point_edge - 
%
%   
%   coords are the barycentric coordinates in the triangle of the new point
%   sum(coords) = 1 and coords(i) > 0
%
%   Copyright (c) 2010 Fethallah Benmansour

if(sum(coords(:)) ~= 1 || min(coords(:)) <= 0)
    error('barycentric coordinates should be strictly in the edge')
end

point_to_add = coords(1)*vertex(:, edge(1)) + ...
               coords(2)*vertex(:, edge(2));

vertex(:, end+1) = point_to_add;

% find if edge blongs to two triangles ?
neigh1 = connectivity(edge(1)).neighbours_idx;
neigh2 = connectivity(edge(2)).neighbours_idx;
k = 0;

opposite_vertices_edge = [0 0];
for i=1:length(neigh1)
    for j=1:length(neigh2)
        if neigh1(i) == neigh2(j)
            k = k + 1;
            opposite_vertices_edge(k) = neigh2(j);
        end
    end
end
if(k == 0 || k > 2)
    error('problem with mesh connectivity');
end

C_point.nb_neighbors = 2+k;
if k==1
    C_point.neighbours_idx = [edge(1) opposite_vertices_edge edge(2)];
else
    C_point.neighbours_idx = [edge(1) opposite_vertices_edge(1)...
                              edge(2) opposite_vertices_edge(2)];
end

connectivity(end+1) = C_point;

modified = [];

for i=1:2

    C_p = connectivity(edge(i));
	Modi= C_p;
    Modi.point_index = edge(i);
    modified = [modified Modi];
    opposite_vertex = edge(mod(i,2)+1);% opposite vertex on edge
    Tab = C_p.neighbours_idx;
    Tab(Tab == opposite_vertex) = length(vertex);
    C_p.neighbours_idx = Tab;
    connectivity(edge(i)) = C_p;
end

for i=1:length(opposite_vertices_edge)
    C_p = connectivity(opposite_vertices_edge(i));
    Modi= C_p;
    Modi.point_index = opposite_vertices_edge(i);
    modified = [modified Modi];
    Tab = insert_point(length(vertex), edge,...
        C_p.neighbours_idx);
    C_p.nb_neighbors = C_p.nb_neighbors + 1;
    C_p.neighbours_idx = Tab;
    connectivity(opposite_vertices_edge(i)) = C_p;
end


return;

% insert_point - 
function Tab = insert_point(point, opposite_edge, neighbords)
    
k = find(neighbords == opposite_edge(1));
tran_to_end = neighbords;
tran_to_end(end-k+1:end) = neighbords(1:k);
tran_to_end(1:end-k)   = neighbords(k+1:end);
Tab = [tran_to_end 0];
if    (tran_to_end(1)     == opposite_edge(2) )
    Tab(end)   = point;
elseif(tran_to_end(end-1) == opposite_edge(2))
    Tab(end)   = Tab(end-1);
    Tab(end-1) = point;
else
    error('should not happen');
end
return;