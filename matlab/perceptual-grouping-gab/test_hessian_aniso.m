% test for anisotropic triangulation with hessian

name = 'lena';
n = 256;
M = load_image(name);
M = rescale( crop(M,n) );

options.order = 2;
H = compute_hessian(M, options);

% one has |l1|>|l2|
[e1,e2,l1,l2] = perform_tensor_decomp(H, 'abs');

% smooth a little
T = perform_tensor_recomp(e1,e2, abs(l2), abs(l1) );
T = perform_blurring(T,2);


plot_tensor_field(T,M);

U = perform_tensor_mapping(T,+1);

% modify here energy and anisotropy


% aspect ratio
[e1,e2,l1,l2] = perform_tensor_decomp(T, 'abs');
alpha = sqrt(abs(l1)./abs(l2));

T = perform_tensor_decomp(e1,e2,1./sqrt(abs(l1)),1./sqrt(abs(l2)));