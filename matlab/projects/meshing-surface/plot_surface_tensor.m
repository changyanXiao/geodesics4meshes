function plot_surface_tensor(vertex,faces,U,landmarks, options)

if not(iscell(U))
    U = {U};
end

lw = getoptions(options, 'lw', 1.5);

col = {'r' 'g' 'b' 'k'};
hold on;
plot_mesh(vertex,faces, options);
for i=1:length(U)
    % plot tensor field
    v = vertex(:,landmarks);
    rho = .05;
    w = v + rho*U{i}(:,landmarks);
    h = plot3( [v(1,:);w(1,:)], [v(2,:);w(2,:)], [v(3,:);w(3,:)], col{i});
    set(h, 'LineWidth', lw);
end