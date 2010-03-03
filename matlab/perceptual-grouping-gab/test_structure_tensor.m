% test for structure tensor

path(path, 'data/');
n = 256;
name = 'curvy-dots';

M = rescale(load_image(name, n));


% compute structure tensor around the curve
sigma1 = 1.5;   % for gradient computation
sigma2 = 30;     % for filtering the tensor field
T = compute_structure_tensor(M,sigma1,sigma2);


[e1,e2,l1,l2] = perform_tensor_decomp(T);
T1 = perform_tensor_recomp(e1,e2,l1*0+1,l2*0+.1);
clf;
plot_tensor_field(T1, M);

U = perform_tensor_mapping(T,+1);