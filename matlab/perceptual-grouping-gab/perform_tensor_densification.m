function T = perform_tensor_densification(T, location, options)

% perform_tensor_densification - interpolate by diffusion tensors
%
% T = perform_tensor_densification(T, location, options);
%
%   interpolate the tensor field by imposing only the value
%   at points where location==1
%
%   You can impose zero tensor in location where options.mask==1
%
%   Copyright (c) 2008 Gabriel Peyre

n = size(T,1);

mask = getoptions(options, 'mask', ones(n));

T0 = T;

niter = getoptions(options, 'niter', 500);

slist = linspace(10,3,niter);

for i=1:niter
    progressbar(i,niter);
    % blurring
    sigma = slist(i);
    T = perform_blurring(T, sigma);
    % impose known values near the curves
    for s=1:3
        T0s = T0(:,:,s); Ts = T(:,:,s);
        Ts(location==1) = T0s(location==1);
        Ts(mask==0) = 0;
        T(:,:,s) = Ts;
    end
    % impose lots of anisotropy
%    U = perform_tensor_mapping(T,+1);
%    h = linspace(0,1,n^2).^2;
%    U(:,:,2) = perform_histogram_equalization(U(:,:,2), h);
%    T = perform_tensor_mapping(U,-1);
end