% test for building an anisotropic metric.

add_base_paths;

rep = '../../results/domain-anisotropy/';
if not(exist(rep))
    mkdir(rep);
end

%%
% load mask

n = 300;
name = 'sala';
M = rescale(load_image(name,n));
M = perform_blurring(M,5)>.5;
M = conv2(double(M), ones(3), 'same')>0;
M = double(M);

%%
% Building tensor field.

options.niter = 500;
H = compute_tensor_domain(M,options);

%%
% Enforce anisotropy

[e1,e2,l1,l2] = perform_tensor_decomp(H);
a = l1+l2;
% l1 = l1./a; l2 = l2./a;
l1 = l1*0+5;
l2 = l2*0+1;
l1(M==0) = 0;
l2(M==0) = 0;
H1 = perform_tensor_decomp(e1,e2,l1,l2);

%%
% Compute tensor eigendecomp.

options.sub = round(n/40);
clf;
plot_tensor_field(H1,1-M, options);
saveas(gcf, [rep name '-metric.png'], 'png');
saveas(gcf, [rep name '-metric.eps'], 'epsc');