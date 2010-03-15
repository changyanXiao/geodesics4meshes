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
H = compute_tensor_domain(M,'anisotropic',options);

%%
% Perform sampling.


%%
% Display.

if strcmp(metric_type, 'anisotropic')
    options.sub = round(n/40);
    clf;
    plot_tensor_field(H,1-M, options);
    saveas(gcf, [rep name '-metric.png'], 'png');
end
saveas(gcf, [rep name '-metric.eps'], 'epsc');