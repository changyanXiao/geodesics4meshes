% test for poisson resolution
% also trace the linking points->boundary provided by the gradient descent
% of the solution

path(path, 'poisson/');

rep  = 'data/';
name = 'toto';
name = 'mm';

repimg = 'results/';
if not(exist(repimg))
    mkdir(repimg);
end

n = 300;
Ma = load_image([rep name],n-10);
M = zeros(n,n,3);
M(6:n-5,6:n-5,:) = Ma;

M = sum(M,3)/3;
mask = 1-(M==M(1));

%% compute the skeleton
[skg,rad] = skeleton(M);
Sk = skg>20;


U = GMGmain(mask);

G = compute_grad(U);
nG = sqrt(sum(G.^2,3));
nG(nG<1e-6) = 1e-6;
G = G./repmat(nG,[1 1 2]);
% divergence, for corners
g1 = compute_grad(G(:,:,1));
g2 = compute_grad(G(:,:,2));
Psi = -g1(:,:,1)-g2(:,:,2);
% rectification for skeleton
Psi2 = U./nG.*Psi;

clf; 
for i=1:4
    subplot(2,2,i);
    imagesc(Psi2>i);
    axis image; axis off;
end
colormap gray(256);
saveas(gcf, [repimg name '-poisson-skeletons'], 'png');

if 0
    % trying to figure out a good rescaling
G = compute_grad(U);
nG = sqrt(sum(G.^2,3));
nG(nG<1e-6) = 1e-6;
G = G./repmat(nG,[1 1 2]);
nG = clamp(nG/std(nG(:)), 0.2, 2);
G = G.*repmat(nG,[1 1 2]);
end

% perform a gradient descent for the point inside
I = find(mask);
[x,y] = ind2sub(size(mask),I);
points = [x(:)';y(:)'];
sel = randperm(size(points,2)); sel = sel(1:1000);
points0 = points(:,sel);

options.dt = 0.2;
options.niter = 200;
traject = compute_vf_trajectory(points0, -G, options);


clf;
subplot(1,2,1);
imagesc(U');
axis image; axis off
title('poisson potential');
colormap jet(256);

subplot(1,2,2);
plot_vf(G(1:4:end,1:4:end,:), mask');
axis tight; axis ij; axis off;
title('vector field');
saveas(gcf, [repimg name '-poisson-vectorfield.png'], 'png');

% plot the trajectories
clf;
hold on;
imagesc(mask'); axis image; axis off; axis ij;
plot_scattered(points0);
plot( squeeze(traject(1,:,:)), squeeze(traject(2,:,:)) );
hold off;
title('association field');
saveas(gcf, [repimg name '-poisson-associations.png'], 'png');