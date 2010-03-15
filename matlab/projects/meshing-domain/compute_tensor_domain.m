function H = compute_tensor_domain(M,metric_type,options)

% compute_tensor_domain - compute a tensor field that follows the normal
%
%   H = compute_tensor_domain(M, metric_type, options);
%
%   M is a binary image, M==1 is the shape.
%   H is a tensor that is rank-1 following the tangent to the domain, and
%   smoothly interpolated inside the shape (and also outside).
%
%   options.niter controls iterations for diffusion (lots are needed !)
%   options.verb controls display of iterations.
%
%   Copyright (c) 2010 Gabriel Peyre

options.null = 0;

n = size(M,1);

switch metric_type
    case 'constant'
        H = zeros(n,n,2,2);
        H(:,:,1,1) = 1; H(:,:,2,2) = 1;
    case 'isotropic'
        H = compute_tensor_domain_iso(M,metric_type,options);
    case 'anisotropic'
        H = compute_tensor_domain_aniso(M,metric_type,options);
    otherwise 
        error('Unknown metric');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = compute_tensor_domain_iso(M,metric_type,options)

options.null = 0;
aniso = getoptions(options, 'isoratio', 5);
niter = getoptions(options, 'niter', 1000);
verb = getoptions(options, 'verb', 1);
n = size(M,1);


B = find(conv2(M, ones(3)/9, 'same')<.999 & M==1);

%%
% compute curvature

T0 = M;
T0 = perform_blurring(T0, 10);
T0 = perform_blurring(T0, 10);
T0 = perform_blurring(T0, 10);

T = zeros(n) + mean(T0(B));
sel1 = [2:n n]; sel2 = [1 1:n-1];
for i=1:niter
    if verb
        progressbar(i,niter);    
    end
    T = ( T+T(sel1,:,:)+T(sel2,:,:)+T(:,sel1,:)+T(:,sel2,:) )/5;
    T(B) = T0(B);
end

T = rescale(T, 1, isoratio);
H = zeros(n,n,2,2);
H(:,:,1,1) = t.^2;
H(:,:,2,2) = t.^2;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = compute_tensor_domain_aniso(M,metric_type,options)

options.null = 0;
aniso = getoptions(options, 'aniso', 5);
niter = getoptions(options, 'niter', 1000);
verb = getoptions(options, 'verb', 1);
n = size(M,1);

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

T = T0;
sel1 = [2:n n]; sel2 = [1 1:n-1];
B = [b; b+n^2; b+2*n^2];
for i=1:niter
    if verb
        progressbar(i,niter);    
    end
    T = ( T+T(sel1,:,:)+T(sel2,:,:)+T(:,sel1,:)+T(:,sel2,:) )/5;
    T(B) = T0(B);
end
H = zeros(n,n,2,2);
H(:,:,1,1) = T(:,:,1);
H(:,:,2,2) = T(:,:,2);
H(:,:,1,2) = T(:,:,3);
H(:,:,2,1) = T(:,:,3);

%%
% Enforce anisotropy

[e1,e2,l1,l2] = perform_tensor_decomp(H);
a = l1+l2;
% l1 = l1./a; l2 = l2./a;
l1 = l1*0+1;
l2 = l2*0+aniso;
l1(M==0) = 0;
l2(M==0) = 0;
H = perform_tensor_decomp(e1,e2,l1,l2);