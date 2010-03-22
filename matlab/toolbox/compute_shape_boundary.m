function bound = compute_shape_boundary(M, trimbound)

% compute_shape_boundary - extract boundary points of a shape
%
% bound = compute_boundary(M, trimbound);
%
%   If trimbound==1, bound is the boundary of the largest connected component of the shape
%   represented by M>mean(M(:))
%
%   Copyright (c) 2009 Gabriel Peyre

if nargin<2
    trimbound = 1;
end

M = double(M>mean(M(:)));
c = contourc(M,[.5 .5]);

b = {}; bsize = [];
while not(isempty(c))
    bsize(end+1) = c(2,1);
    b{end+1} = c(:,2:bsize(end)+1);
    c(:,1:bsize(end)+1) = [];
end
    
if trimbound
    [tmp,I] = max(bsize);
    bound = b{I};
else
    bound = [b{:}];
end

bound = bound(2:-1:1,:);