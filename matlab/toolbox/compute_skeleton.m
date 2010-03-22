function [SK,DT,LFS] = compute_skeleton(M, tau)

% compute_skeleton - compute skeleton
%
%   [SK,DT,LFS] = compute_skeleton(M,tau)
%
%   SK is the skeleton map.
%   DT is the distance transform of the boundary of the shape.
%   LFS is the local feature size, the distance to the skeletton.
%
%   tau is a small threshold for the skeleton.
%
%   Copyright (c) 2010 Gabriel Peyre

n = size(M,1);


if 1
    pstart = compute_shape_boundary(M, 0);
else
    B = conv2(M, ones(3)/9, 'same')<.999 & M==1;
    [x,y] = ind2sub(size(M), find(B));
    pstart = [x(:)'; y(:)'];
end

[DT,S,Q] = perform_fast_marching(ones(n), pstart);

nbound = size(pstart,2);
G = grad(Q);
G(G<-nbound/2) = G(G<-nbound/2) + nbound;
G(G>nbound/2) = G(G>nbound/2) - nbound;
G = sqrt(sum(G.^2,3));

if nargin<2
    tau = nbound/30;
end


SK = G>tau;
SK(M==0) = 0;

if nargout>2
    [x,y] = ind2sub(size(M), find(SK));
    pstart = [x(:)'; y(:)'];
    [LFS,S,Q] = perform_fast_marching(ones(n), pstart);
end
