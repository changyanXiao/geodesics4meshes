function le = compute_tensor_distance(v, t)

%% compute_tensor_distance - compute the L2 metric for the tensor
%
%   le = compute_tensor_distance(v, t);
%
%   le^2 = v'*t*v
%
%   t can be in 3x3 format or in 6x1 format.
%
%   Copyright (c) Gabriel Peyre


if size(t,1)==6
    m = size(t,2);
    t = reshape(t, [6 1 m]);
    t1 = zeros(3,3,m);
    t1(1,1,:) = t(1,1,:);
    t1(2,2,:) = t(2,1,:);
    t1(3,3,:) = t(3,1,:);
    t1(1,2,:) = t(4,1,:);
    t1(2,3,:) = t(5,1,:);
    t1(1,3,:) = t(6,1,:);
    t1(2,1,:) = t1(1,2,:);
    t1(3,1,:) = t1(1,3,:);
    t1(3,2,:) = t1(2,3,:);
    t = t1;
end

m = size(t,3);
le = zeros(m,1);
for i=1:m
    le(i) = v(:,i)'*t(:,:,i)*v(:,i);
end
le = sqrt(le);