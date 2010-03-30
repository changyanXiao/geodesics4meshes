function [vertex, connectivity, modified] = add_point_triangle(vertex, connectivity, triangle, coords)

% add_point_triangle - 
%
%   
%   coords are the barycentric coordinates in the triangle of the new point
%   sum(coords) = 1 and coords(i) > 0
%
%   Copyright (c) 2010 Fethallah Benmansour

if(sum(coords(:)) ~= 1 || min(coords(:)) <= 0)
    error('barycentric coordinates should be strictly in the triangle')
end

point_to_add = coords(1)*vertex(:, triangle(1)) + ...
               coords(2)*vertex(:, triangle(2)) + ...
               coords(3)*vertex(:, triangle(3)) ;

           
vertex(:, end+1) = point_to_add;

C_point.nb_neighbors = 3;
C_point.neighbours_idx = triangle;

connectivity(end+1) = C_point;
modified = [];
for i=1:3
    C_p = connectivity(triangle(i));
    opposite_edge = [triangle(mod(i,3)+1), triangle(mod(i+1,3)+1)];
    Modi= C_p;
    Modi.point_index = triangle(i);
    modified = [modified Modi];
	C_p.nb_neighbors = C_p.nb_neighbors + 1;
    neighbords = C_p.neighbours_idx;
	TAB = insert_point(length(vertex), opposite_edge, neighbords);
    C_p.neighbours_idx = TAB;
    connectivity(triangle(i)) = C_p;
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
