% test for anisotropic propagation

path(path, 'toolbox/');

n = 128;

% load smooth vector field
v = perform_vf_normalization(perform_blurring(randn(n,n,2), 50));
w = cat(3,v(:,:,2),-v(:,:,1));

pstart = [n/2;n/2];
Dlist = [10 20 30 50 80 1e9];

lambda = .05;
T = perform_tensor_decomp(v,w,ones(n), ones(n)*lambda);
    
clf;
for i=1:length(Dlist);
    Dmax = Dlist(i); 
    [U,dUx,dUy,V,L] = fm2dAniso([1;1], T, pstart, Dmax);    
    subplot(2,3,i);
    V = U; V(V>1e8) = 0;
    imageplot(V);
    colormap jet(256);
end