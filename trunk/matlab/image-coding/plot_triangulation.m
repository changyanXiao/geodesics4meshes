function plot_triangulation(vertex,faces, M, options)

% plot_triangulation - plot a 2D triangulation
%
%   plot_triangulation(vertex,face,M, options);
%
%   The point are assume to be in (1,n).
%
%   Copyright (c) 2009 Gabriel Peyre

options.null = 0;
lw = getoptions(options, 'edgewidth', 2);
ms = getoptions(options, 'vertexsize', 20);

edges = compute_edges(faces);

vertex = vertex(2:-1:1,:);

if not(isempty(M))
    imageplot(M);
end
hold on;
hh = plot_edges(edges, vertex, 'k');
set(hh, 'LineWidth', lw);
hh = plot(vertex(1,:),vertex(2,:), 'k.');
set(hh, 'MarkerSize', ms);
hold off;