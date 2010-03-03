% test for dense tensor field interpolation

path(path, 'data/');
name = 'two-lines';
name = 'sparse-curve';

rep = 'results/tensor-dense/';
if not(exist(rep))
    mkdir(rep);
end

n = 200;
M = load_image(name, n);
M = rescale(clamp(M,0,255));

% compute structure tensor around the curve
sigma1 = 1.5;   % for gradient computation
sigma2 = 6;     % for filtering the tensor field
T0 = compute_structure_tensor(M,sigma1,sigma2);


T = T0;
U = perform_tensor_mapping(T,+1);
Theta = U(:,:,3);
s = 1e-5;
U(:,:,1) = U(:,:,1) .* (U(:,:,1)>s);
U(:,:,2) = .5;
T1 = perform_tensor_mapping(U,-1);

% display tensor field
clf;
plot_tensor_field(T1, M);
saveas(gcf, [rep name '-sparse-tf.png'], 'png');


% compute eigenvectors
[e1,e2,l1,l2] = perform_tensor_decomp(T);

% set to 0/1 the eigenvalues
T = perform_tensor_recomp(e1,e2,l1*0+1,l2*0);
% to avoid using 2x2 matrices
T = cat(3, T(:,:,1,1), T(:,:,2,2), T(:,:,1,2));
T0 = T;

% perform iterative diffusion
A = M<=.01;

[Y,X] = meshgrid(1:n,1:n);
n1 = 65;
n2 = 155;
if strcmp(name, 'two-lines')
    A = A.*(abs(Theta-Theta(85,30))<.01) + ...
        A.*(abs(Theta-Theta(105,174))<.01);
end

niter = 500;
slist = linspace(15,3,niter);
for i=1:niter
    progressbar(i,niter);
    % blurring
    sigma = slist(i);
    T = perform_blurring(T, sigma);
    % impose known values near the curves
    for s=1:3
        T0s = T0(:,:,s); Ts = T(:,:,s);
        Ts(A==1) = T0s(A==1);
        T(:,:,s) = Ts;
    end
    % impose lots of anisotropy
%    U = perform_tensor_mapping(T,+1);
%    h = linspace(0,1,n^2).^2;
%    U(:,:,2) = perform_histogram_equalization(U(:,:,2), h);
%    T = perform_tensor_mapping(U,-1);
end


U = perform_tensor_mapping(T,+1);
U(:,:,2) = .5;
T1 = perform_tensor_mapping(U,-1);


clf;
plot_tensor_field(T1, M);
saveas(gcf, [rep name '-dense-tf-' num2str(k) '.png'], 'png')