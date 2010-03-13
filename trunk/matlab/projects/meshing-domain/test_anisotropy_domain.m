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
% Compute normals

G0 = grad(perform_blurring(M,5));

%% 
% Tensor it.

T0 = [];
T0(:,:,1) = G0(:,:,2).^2;
T0(:,:,2) = G0(:,:,1).^2;
T0(:,:,3) = -G0(:,:,1).*G0(:,:,2);

%%
% Boundary 

b = find(conv2(M, ones(3)/9, 'same')<.999 & M==1);

%%
% Diffuse it.

niter = 1000;
T = T0;
sel1 = [2:n n]; sel2 = [1 1:n-1];
B = [b; b+n^2; b+2*n^2];
for i=1:niter
    progressbar(i,niter);    
    T = ( T+T(sel1,:,:)+T(sel2,:,:)+T(:,sel1,:)+T(:,sel2,:) )/5;
    T(B) = T0(B);
end

%%
% Compute tensor eigendecomp.

H = zeros(n,n,2,2);
H(:,:,1,1) = T(:,:,1);
H(:,:,2,2) = T(:,:,2);
H(:,:,1,2) = T(:,:,3);
H(:,:,2,1) = T(:,:,3);

%%
% Display.

[e1,e2,l1,l2] = perform_tensor_decomp(H);
a = l1+l2;
% l1 = l1./a; l2 = l2./a;
l1 = l1*0+5;
l2 = l2*0+1;
l1(M==0) = 0;
l2(M==0) = 0;
H1 = perform_tensor_decomp(e1,e2,l1,l2);
options.sub = round(n/40);
clf;
plot_tensor_field(H1,1-M, options);
saveas(gcf, [rep name '-metric.png'], 'png');
saveas(gcf, [rep name '-metric.eps'], 'epsc');

return;

s = 2;
clf;
plot_vf(e1(1:s:end,1:s:end, :), M);