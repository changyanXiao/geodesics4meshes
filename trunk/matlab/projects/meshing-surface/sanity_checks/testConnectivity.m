%%
% Test connectivity to matlab.
clear all; close all; clc;

add_base_paths;
%%
name = 'cylinder';
options.name = name;
[vertex,faces] = read_mesh(name);
%%
tic
Connectivity = ComputeMeshConnectivity(vertex-1, faces-1);
toc
%% sanity check for the connectivity
figure;
hold on;
point_idx = 300;
colors = jet(Connectivity(point_idx).nb_neighbors);
plot3(vertex(2, point_idx), vertex(1, point_idx), vertex(3, point_idx), 'or', 'MarkerSize', 5);
plot3(vertex(2, Connectivity(point_idx).neighbours_idx(1:end)),...
	  vertex(1, Connectivity(point_idx).neighbours_idx(1:end)),...
      vertex(3, Connectivity(point_idx).neighbours_idx(1:end)), '-k', 'linewidth', 2);
plot3([vertex(2, Connectivity(point_idx).neighbours_idx(1)), vertex(2, Connectivity(point_idx).neighbours_idx(end))],...
      [vertex(1, Connectivity(point_idx).neighbours_idx(1)), vertex(1, Connectivity(point_idx).neighbours_idx(end))],...
      [vertex(3, Connectivity(point_idx).neighbours_idx(1)), vertex(3, Connectivity(point_idx).neighbours_idx(end))], '-k');  
for k=1:Connectivity(point_idx).nb_neighbors
    plot3([vertex(2, point_idx) vertex(2, Connectivity(point_idx).neighbours_idx(k))],...
          [vertex(1, point_idx) vertex(1, Connectivity(point_idx).neighbours_idx(k))],...
          [vertex(3, point_idx) vertex(3, Connectivity(point_idx).neighbours_idx(k))], '-k', 'lineWidth', 3, 'color', colors(k,:) );
end
plot_mesh(vertex, faces);
%% sanity check for adding a point in triangle
triangle = faces(:, 300);
coords = [0.5 0.3 0.2];
point_to_add = coords(1)*vertex(:, triangle(1)) + ...
               coords(2)*vertex(:, triangle(2)) + ...
               coords(3)*vertex(:, triangle(3)) ;
[ver, connec, mod] = add_point_triangle(vertex, Connectivity, triangle, coords);
%%
figure;
hold on;
axis off,
set(gca, 'color', 'none')
for i=1:3
    point_idx = triangle(i);
    plot3(ver(2, point_idx), ver(1, point_idx), ver(3, point_idx),  'or', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot3(ver(2, connec(point_idx).neighbours_idx(1:end)),...
        ver(1, connec(point_idx).neighbours_idx(1:end)),...
        ver(3, connec(point_idx).neighbours_idx(1:end)), '-k', 'linewidth', 2);
    plot3([ver(2, connec(point_idx).neighbours_idx(1)), ver(2, connec(point_idx).neighbours_idx(end))],...
        [ver(1, connec(point_idx).neighbours_idx(1)), ver(1, connec(point_idx).neighbours_idx(end))],...
        [ver(3, connec(point_idx).neighbours_idx(1)), ver(3, connec(point_idx).neighbours_idx(end))], '-k');
    for k=1:connec(point_idx).nb_neighbors
        plot3([ver(2, point_idx) ver(2, connec(point_idx).neighbours_idx(k))],...
            [ver(1, point_idx) ver(1, connec(point_idx).neighbours_idx(k))],...
            [ver(3, point_idx) ver(3, connec(point_idx).neighbours_idx(k))], '-', 'lineWidth', 3, 'color', 'k' );
    end
end
point_idx = length(ver);% added vertex is at the end 

colors = jet(connec(point_idx).nb_neighbors);
plot3(ver(2, point_idx), ver(1, point_idx), ver(3, point_idx), 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot3(ver(2, connec(point_idx).neighbours_idx(1:end)),...
	  ver(1, connec(point_idx).neighbours_idx(1:end)),...
      ver(3, connec(point_idx).neighbours_idx(1:end)), '-k', 'linewidth', 2);
plot3([ver(2, connec(point_idx).neighbours_idx(1)), ver(2, connec(point_idx).neighbours_idx(end))],...
      [ver(1, connec(point_idx).neighbours_idx(1)), ver(1, connec(point_idx).neighbours_idx(end))],...
      [ver(3, connec(point_idx).neighbours_idx(1)), ver(3, connec(point_idx).neighbours_idx(end))], '-k');  
for k=1:connec(point_idx).nb_neighbors
    plot3([ver(2, point_idx) ver(2, connec(point_idx).neighbours_idx(k))],...
          [ver(1, point_idx) ver(1, connec(point_idx).neighbours_idx(k))],...
          [ver(3, point_idx) ver(3, connec(point_idx).neighbours_idx(k))], '-k', 'lineWidth', 3);
end
view(-150, 20);
axis off,
set(gca, 'color', 'none')
%% now check a neighbor
figure;
hold on;
point_idx = triangle(i);
plot3(ver(2, point_idx), ver(1, point_idx), ver(3, point_idx),  'or', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(ver(2, Connectivity(point_idx).neighbours_idx(1:end)),...
    ver(1, Connectivity(point_idx).neighbours_idx(1:end)),...
    ver(3, Connectivity(point_idx).neighbours_idx(1:end)), '-k', 'linewidth', 2);
plot3([ver(2, Connectivity(point_idx).neighbours_idx(1)), ver(2, Connectivity(point_idx).neighbours_idx(end))],...
    [ver(1, Connectivity(point_idx).neighbours_idx(1)), ver(1, Connectivity(point_idx).neighbours_idx(end))],...
    [ver(3, Connectivity(point_idx).neighbours_idx(1)), ver(3, Connectivity(point_idx).neighbours_idx(end))], '-k');
for k=1:Connectivity(point_idx).nb_neighbors
    plot3([ver(2, point_idx) ver(2, Connectivity(point_idx).neighbours_idx(k))],...
        [ver(1, point_idx) ver(1, Connectivity(point_idx).neighbours_idx(k))],...
        [ver(3, point_idx) ver(3, Connectivity(point_idx).neighbours_idx(k))], '-', 'lineWidth', 3, 'color', 'k' );
end
view(-150, 20);
axis off,
set(gca, 'color', 'none')
title('before')

figure;
hold on;
point_idx = triangle(i);
plot3(ver(2, point_idx), ver(1, point_idx), ver(3, point_idx),  'or', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(ver(2, connec(point_idx).neighbours_idx(1:end)),...
    ver(1, connec(point_idx).neighbours_idx(1:end)),...
    ver(3, connec(point_idx).neighbours_idx(1:end)), '-k', 'linewidth', 2);
plot3([ver(2, connec(point_idx).neighbours_idx(1)), ver(2, connec(point_idx).neighbours_idx(end))],...
    [ver(1, connec(point_idx).neighbours_idx(1)), ver(1, connec(point_idx).neighbours_idx(end))],...
    [ver(3, connec(point_idx).neighbours_idx(1)), ver(3, connec(point_idx).neighbours_idx(end))], '-k');
for k=1:connec(point_idx).nb_neighbors
    plot3([ver(2, point_idx) ver(2, connec(point_idx).neighbours_idx(k))],...
        [ver(1, point_idx) ver(1, connec(point_idx).neighbours_idx(k))],...
        [ver(3, point_idx) ver(3, connec(point_idx).neighbours_idx(k))], '-', 'lineWidth', 3, 'color', 'k' );
end
view(-150, 20);
axis off,
set(gca, 'color', 'none')
title('after')
%% Sanity check add point to edge
edge = triangle(1:2);
coords = [0.5 0.5];
[ver, connec, mod] = add_point_edge(vertex, Connectivity, edge, coords);
%%
figure;
hold on;
axis off,
set(gca, 'color', 'none')
for i=1:2
point_idx = edge(i);
plot3(ver(2, point_idx), ver(1, point_idx), ver(3, point_idx),  'or', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(ver(2, connec(point_idx).neighbours_idx(1:end)),...
	  ver(1, connec(point_idx).neighbours_idx(1:end)),...
      ver(3, connec(point_idx).neighbours_idx(1:end)), '-k', 'linewidth', 2);
plot3([ver(2, connec(point_idx).neighbours_idx(1)), ver(2, connec(point_idx).neighbours_idx(end))],...
      [ver(1, connec(point_idx).neighbours_idx(1)), ver(1, connec(point_idx).neighbours_idx(end))],...
      [ver(3, connec(point_idx).neighbours_idx(1)), ver(3, connec(point_idx).neighbours_idx(end))], '-k');  
for k=1:connec(point_idx).nb_neighbors
    plot3([ver(2, point_idx) ver(2, connec(point_idx).neighbours_idx(k))],...
          [ver(1, point_idx) ver(1, connec(point_idx).neighbours_idx(k))],...
          [ver(3, point_idx) ver(3, connec(point_idx).neighbours_idx(k))], '-', 'lineWidth', 3, 'color', 'k' );
end
point_idx = length(ver);% added vertex is at the end 
colors = jet(connec(point_idx).nb_neighbors);
plot3(ver(2, point_idx), ver(1, point_idx), ver(3, point_idx), 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
view(-180, 20);
end