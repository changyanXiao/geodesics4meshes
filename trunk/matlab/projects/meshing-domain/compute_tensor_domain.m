function H = compute_tensor_domain(M,options)

% compute_tensor_domain - compute a tensor field that follows the normal
%
%   H = compute_tensor_domain(M, options);
%
%   M is a binary image, M==1 is the shape.
%   H is a tensor that is rank-1 following the tangent to the domain, and
%   smoothly interpolated inside the shape (and also outside).
%
%   options.niter controls iterations for diffusion (lots are needed !)
%   options.verb controls display of iterations.
%
%   Copyright (c) 2010 Gabriel Peyre

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