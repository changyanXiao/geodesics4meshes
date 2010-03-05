% test for meshing the inside of an object

rep  = 'data/';
name = 'mm';
n = 220;
Ma = load_image([rep name],n-10);
M = zeros(n,n,3);
M(6:n-5,6:n-5,:) = Ma;

repimg = 'results/meshing/';
if ~exist(repimg)
    mkdir(repimg);
end


M1 = sum(M,3);
mask = 1-(M1==M1(1));
p = 1000;

%% compute distance to boundary field
% compute boundary points
h = ones(3)/9;
H = perform_convolution(double(mask),h);
B = H>2/9 & H<1-2/9;
I = find(B);
[x,y] = ind2sub(size(M),I);
boundary = [x(:)';y(:)'];
% compute distance to boundary
[D,Z,Q] = perform_fast_marching(ones(n), boundary);


use_adaptive = 1;
if use_adaptive
    R = 0.8;
    D1 = min(rescale(D),R);
    H = sqrt( R^2 - (D1-R).^2 ) * n;
    W = rescale( D, 0.1,1 );
else
    W = ones(n);
end

if 0
    %% perform sampling using farthest point
    L = zeros(n) - Inf;
    I = find(mask); L(I) = Inf;
    vertex = [n/2;n/2];
    options.constraint_map = L;
    vertex = farthest_point_sampling(W, vertex, p, options );

    %% compute the associated triangulation
    [D,Z,Q] = perform_fast_marching(W, vertex, options);
    faces = compute_voronoi_triangulation(Q, vertex);

    %% display
    clf;
    hold on;
    imagesc(rescale(M)); axis image; axis off;
    plot(vertex(2,:), vertex(1,:), 'b.', 'MarkerSize', 20);
    plot_edges(compute_edges(faces), vertex(2:-1:1,:), 'r');
    hold off;
    axis tight; axis image; axis off;
    colormap gray(256);
    axis ij;

    str = [name '-mesh-' num2str(p)];
    if use_adaptive
        str = [str '-adaptive'];
    end
    saveas(gcf, [repimg str '.png'], 'png');
end


%% compute a smoothed boundary curve
[D,Z,Q] = perform_fast_marching(ones(n), boundary);
% signed distance
D1 = D .* (mask-0.5);
c = perform_contour_extraction(D1, 0);
% smooth the curve
h = ones(31,1); h = h/sum(h(:));
for i=1:2
    c(i,:) = perform_convolution(c(i,:), h);
end
% turn into unique point
boundary1 = round( c*(n-1)+1 );
[tmp,I,J] = unique(boundary1', 'rows');
I = sort(I);
boundary1 = boundary1(:,I);

%% compute a curvilinear 2D indexing using (curv.abs,distance)
lambda = 1;
boundary1 = rescale(boundary1);
s = compute_cuvilinear_abscice( boundary1 );
[DD,Q] = perform_distance_transform(ones(n*lambda), boundary1);

% [DD,Z,Q] = perform_fast_marching(ones(n), boundary1);
s = rescale(s(Q));                  % absice
t = rescale(H); t(mask==0)=-Inf;    % distance
% number of sampling points along each dimension
ns = 200; nt = 10;
[T,S] = meshgrid(linspace(0,1,nt), linspace(0,1,ns));
% inverse the mapping to find the sampling locations
vert = zeros(2,ns*nt);
for i=1:ns
    for j=1:nt
        progressbar( (i-1)*nt+j, (ns*nt) );
        s1 = S(i,j); 
        t1 = T(i,j);
        d = (s-s1).^2 + (t-t1).^2;
        [tmp,I] = min(d(:));
        [x,y] = ind2sub(size(H), I);
        vert(:,i+(j-1)*ns) = [x;y];
    end
end

return;


%% compute the adaptive R
q = max(Q(:));
r = zeros(q,1);
for i = 1:q
    v = D(Q==i & mask==1);
    if not(isempty(v))
        r(i) = max(v);
    end
end
RQ = r(Q);


%% compute height field H
[D,Z,Qtmp] = perform_fast_marching(ones(n), boundary);
% R = 0.2 * D(end/2,end/2); % rescale(D)
D1 = min(D,RQ);
D1 = max(D1,0);
H = sqrt( RQ.^2 - (D1-RQ).^2 )*n;
%% smooth it a bit
Z = H.*mask;
h = ones(21); h = h / sum(h(:));
Z = perform_convolution(Z,h);



%% compute the 3D mesh
% compute the elevation of the vertices
H = Z;
p = size(vertex,2);
vertex3d = zeros(3,p); 
vertex3d(1:2,:) = vertex;
for i=1:p
    vertex3d(3,i) = H(round(vertex(1,i)),round(vertex(2,i)));
end
% symmetrize
vertex3d = [vertex3d, vertex3d];
vertex3d(3,end/2+1:end) = -vertex3d(3,end/2+1:end);
faces3d = [faces; faces+p];
% display (should reorient the faces)
figure;
clf;
plot_mesh(vertex3d,faces3d);
lighting flat;