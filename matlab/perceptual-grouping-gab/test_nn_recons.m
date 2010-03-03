path(path, 'data/');
name = 'sparse-curve';

rep = 'results/tensor-dense/';
if not(exist(rep))
    mkdir(rep);
end

n = 150;
M = load_image(name, n);
M = rescale(clamp(M,0,255));

% compute structure tensor around the curve
sigma1 = 1.5;   % for gradient computation
sigma2 = 6;     % for filtering the tensor field
T1 = compute_structure_tensor(M,sigma1,sigma2);

% compute eigenvectors
[e1,e2,l1,l2] = perform_tensor_decomp(T1);

% display tensor field
%clf;
%plot_tensor_field(T1, M);
%saveas(gcf, [rep name '-sparse-tf.png'], 'png');

% set to 0/1 the eigenvalues
T1 = perform_tensor_recomp(e1,e2,l1*0+1,l2*0);
% to avoid using 2x2 matrices
T1 = cat(3, T1(:,:,1,1), T1(:,:,2,2), T1(:,:,1,2));
T0 = T1;

% perform iterative diffusion
sigma = 3;
niter = 50;
for i=1:niter
    progressbar(i,niter);
    % blurring
    T1 = perform_blurring(T1, sigma);
    % impose known values near the curves
    for s=1:3
        T0s = T0(:,:,s); Ts = T1(:,:,s);
        Ts(M==0) = T0s(M==0);
        T1(:,:,s) = Ts;
    end
    % impose lots of anisotropy
    U = perform_tensor_mapping(T1,+1);
    h = linspace(0,1,n^2).^2;
    U(:,:,2) = perform_histogram_equalization(U(:,:,2), h);
    T1 = perform_tensor_mapping(U,-1);
end

U = perform_tensor_mapping(T1,+1);
U(:,:,1) = 1;
U(:,:,2) = clamp(U(:,:,2), .8,1);
T1 = perform_tensor_mapping(U,-1);

T = ones(n,n,2,2);
T(:,:,1,1) = T1(:,:,1);
T(:,:,2,2) = T1(:,:,2);
T(:,:,1,2) = T1(:,:,3);
T(:,:,2,1) = T1(:,:,3);

% picking points
clf;
nstart = 8;
start_points=ones(2,nstart);

imageplot(M);
axis image; axis off;
colormap gray(256);

hold on;
for nbp = 1:nstart
    disp('Pick a point.');
    start_point = ginput(1);
    plot(start_point(1),start_point(2),'rs');
    start_points(1,nbp)=start_point(1);
    start_points(2,nbp)=start_point(2);
end
hold off;
start_points = round(start_points);

D = perform_fast_marching(T, start_points);


%clf;
%imageplot(M);
%axis image;

%clf;
%options.sub = round(n/15);
%options.color = 'k';
plot_tensor_field(T, D);
hold on;
%imageplot(D');
colormap jet(256);
h = plot(start_points(1,:),start_points(2,:), 'b.');    
set(h, 'MarkerSize', 20);
%colormap jet(256);

%saveas(gcf, [rep name '-aniso-' num2str(ianiso) '.png'], 'png');