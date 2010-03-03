load points1;
load distance;

pts = [107; 137];
pend = [52;121];

% distance to end points
n = size(D,1);
[Y,X] = meshgrid(1:n,1:n);
d = sqrt( (X-pend(1)).^2+(Y-pend(2)).^2 );
d = rescale(d, 0, max(D(:)));

options.method = 'continuous';
options.trim_path = 0;
pathc = compute_geodesic(D,pts,options);
options.method = 'discrete';
pathd = compute_geodesic(D,pts,options);


clf;
hold on;
imagesc(D'); axis image; axis off;
h = plot(pathd(1,:), pathd(2,:), 'r');
set(h, 'LineWidth', 2);
axis image; axis off;
h = plot(pathc(1,:), pathc(2,:), 'k');
set(h, 'LineWidth', 2);
colormap jet;

return;


b = 1; pathc = [];
while b==1
    clf;
    hold on;
    imagesc(D'); axis image; axis off;
    if not(isempty(pathc))
        h = plot(pathc(1,:), pathc(2,:), 'k');
        set(h, 'LineWidth', 2);
        h = plot(pathd(1,:), pathd(2,:), 'r');
        set(h, 'LineWidth', 2);
        axis image; axis off;
    end
    [x,y,b] = ginput(1);
    pts = round([x;y]);
    pts = [137; 107];
    options.method = 'continuous';
    pathc = compute_geodesic(D,pts,options);
    options.method = 'discrete';
    pathd = compute_geodesic(D,pts,options);
end