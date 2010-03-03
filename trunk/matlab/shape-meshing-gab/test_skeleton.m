% test for distance transform and skeletton extraction

rep  = 'data/';
name = 'toto';

n = 350;
Ma = load_image([rep name],n-10);
M = zeros(n,n,3);
M(6:n-5,6:n-5,:) = Ma;

M = sum(M,3)/3;
mask = 1-(M==M(1));

%% compute the skeleton
[skg,rad] = skeleton(M);
Sk = skg>20;

%% compute the distance to skeletton
I = find(Sk);
[x,y] = ind2sub(size(M),I);
Skpos = [x(:)';y(:)'];
[Dsk,Z,Q] = perform_fast_marching(ones(n), Skpos);

%% compute boundary points
h = ones(3)/9;
H = perform_convolution(double(mask),h);
B = H>2/9 & H<1-2/9;
I = find(B);
[x,y] = ind2sub(size(M),I);
boundary = [x(:)';y(:)'];


%% compute a smoothed boundary curve
disp('--> Smooth boundary');
[Dbound,Z,Q] = perform_fast_marching(ones(n), boundary);
% signed distance
D1 = Dbound .* (mask-0.5);
boundary1 = perform_contour_extraction(D1, 0);
% smooth the curve
h = ones(31,1); h = h/sum(h(:));
for i=1:2
    boundary1(i,:) = perform_convolution(boundary1(i,:), h, 'per');
end
m = 1500;
% resample the curve
boundary1 = perform_curve_resampling(boundary1, m);
boundary1 = round( boundary1*(n-1)+1 );
boundary1 = unique(boundary1', 'rows')';

%% compute nearest point
disp('--> Slow distance transform');
% [D,Q] = perform_distance_transform(mask, boundary1);
[D,Z,Q] = perform_fast_marching(ones(n), boundary1);

%% compute radius
Q(Q==0) = 1;
Indb = boundary1(1,Q) + n*(boundary1(2,Q)-1); 
Indb = reshape(Indb,n,n);
R = Dsk(Indb);

%% compute extrudation
hmax = 100*n/300;
r = 0.1;
r1 = max(r-Dbound, 0);
E = sqrt( r.^2 - r1.^2 );
E1 = E.*mask;

r = R;
r1 = max(r-Dbound, 0);
E = sqrt( r.^2 - r1.^2 );
E2 = E.*mask;

rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end


clf;
subplot(2,3,1);
imagesc(Dbound.*mask); axis image; axis off;
title('Distance to boundary');
subplot(2,3,2);
imagesc(Sk + mask); axis image; axis off;
title('Skeleton');
subplot(2,3,3);
imagesc(Dsk.*mask); axis image; axis off;
title('Distance to sk');
subplot(2,3,4);
contour(Q, 150); axis image; axis off; axis ij;
title('Boundary assocations');
subplot(2,3,5);
imagesc(R); axis image; axis off; axis ij;
title('Local radius');
colormap jet(256);

saveas(gcf, [rep name '-distfunc.png'], 'png');

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


saveas(gcf, [rep name '-extrude.png'], 'png');


return;


% compute the max distance function together with the point where it is reached
Dm = zeros(m,1);
S = zeros(m,1);
for i = 1:m
    I = find(Q==i);
    if not(isempty(I))
        v = D(I);
        [Dm(i), k] = max(v);
        S(i) = I(k);
    end
end

Sk = zeros(n);
Sk(S) = 1;

I = find(Q>0);
DM = zeros(n);
DM(I) = Dm(Q(I));