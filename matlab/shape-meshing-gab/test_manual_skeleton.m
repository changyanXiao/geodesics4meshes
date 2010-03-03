% test for shape extrusion using hand-drawing skeleton

% test for distance transform and skeletton extraction

rep  = 'data/';
name = 'toto';

n = 350;
Ma = load_image([rep name],n-10);
M = zeros(n,n,3);
M(6:n-5,6:n-5,:) = Ma;

M = sum(M,3)/3;
mask = 1-(M==M(1));
sk = pick_curves(mask);

%% compute distance to skeleton
I = find(sk);
[x,y] = ind2sub(size(mask),I);
skpoints = [x(:)';y(:)'];
nsk = length(I); % number of skeleton points
Isk = x(:) + (y(:)-1)*n;
CM = zeros(n) + Inf; CM(mask==0) = -Inf;
options.constraint_map = CM;
[Dsk,Z,Qsk] = perform_fast_marching(ones(n), skpoints, options);
Dsk(Dsk==Inf) = max(Dsk(Dsk~=Inf));

if 0
    % use distance transform for the correspondance
    Qsk(Qsk==0) = 1;
    [Dsk,Qsk] = eucdist2(logical(sk));
elseif 0
    % perform a gradient descent for the point inside
    I = find(mask);
    [x,y] = ind2sub(size(mask),I);
    points = [x(:)';y(:)'];
    sel = randperm(size(points,2)); sel = sel(1:10000);
    sel = 1:size(points,2);
    points0 = points(:,sel);
    % compute gradient of the distance function
    G = perform_vf_normalization(compute_grad(Dsk));
    % compute the trajectory to the boundary
    options.dt = 0.2;
    options.niter = 200;
    traject = compute_vf_trajectory(points0, -G, options);
    
    Isk = round(traject(1,end,:)) + n*( round(traject(2,end,:)) - 1);
    Qsk = ones(n); Qsk(mask==1) = Isk;
    
    % plot the trajectories
    if 0
    sel = randperm(size(points,2)); sel = sel(1:5000);
    clf;
    hold on;
    imagesc(mask'); axis image; axis off; axis ij;
    plot_scattered(points0);
    plot( squeeze(traject(1,:,sel)), squeeze(traject(2,:,sel)) );
    hold off;
    title('association field');
    saveas(gcf, [repimg name '-geoddist-associations.png'], 'png');
    end
end

%% compute distance to boundary
h = ones(3)/9;
H = perform_convolution(double(mask),h);
B = H>2/9 & H<1-2/9;
I = find(B);
[x,y] = ind2sub(size(M),I);
boundary = [x(:)';y(:)'];
Ibound = x(:) + (y(:)-1)*n;
[Dbound1,Z,Qbound1] = perform_fast_marching(ones(n), boundary);
[Dbound,Qbound] = eucdist2(logical(1-mask));

%% compute the local radius by geodesic interpolation
f = Dbound(Isk)';
% subsampling
sel = randperm(nsk); sel = sel(1:100);
options.sigma = 1/n;
options.sigma = 20/n;
options.alpha = 2;
options.method = 'gaussian';
Rsk = perform_geodesic_interpolation(ones(n),skpoints(:,sel),f(sel),options);
% Rsk = Dbound((Qsk));

%% compute the extrusion
hmax = 100*n/300;
r = 10;
r1 = max(r-Dbound, 0);
E = sqrt( r.^2 - r1.^2 );
E1 = E.*mask;

r = Rsk;
r1 = max(r-Dbound, 0);
E = sqrt( r.^2 - r1.^2 );
E2 = E.*mask;

repimg = 'results/extrusion/';
if not(exist(repimg))
    mkdir(repimg);
end


clf;
subplot(2,3,1);
imagesc(Dbound.*mask); axis image; axis off;
title('Distance to boundary');
subplot(2,3,2);
imagesc(sk + mask); axis image; axis off;
title('Skeleton');
subplot(2,3,3);
imagesc(Dsk.*mask); axis image; axis off;
title('Distance to sk');
subplot(2,3,4);
imagesc(Qsk.*mask); axis image; axis off; axis ij;
title('Boundary assocations');
subplot(2,3,5);
imagesc(Rsk.*mask); axis image; axis off; axis ij;
title('Local radius');
colormap jet(256);

saveas(gcf, [repimg name '-distfunc.png'], 'png');

figure;
clf;
subplot(1,2,1);
surf(rescale(E1,0,hmax), E*0); 
shading interp;
colormap gray(256);
axis equal; axis off;
view(-137, 60);
camlight;

subplot(1,2,2);
hmax = 100*n/300;
surf(rescale(E2,0,hmax), E*0); 
shading interp;
colormap gray(256);
axis equal; axis off;
view(-137, 60);
camlight;


saveas(gcf, [repimg name '-extrude.png'], 'png');
