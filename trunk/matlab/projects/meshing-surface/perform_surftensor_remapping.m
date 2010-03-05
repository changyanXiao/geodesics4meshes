function T = perform_surftensor_remapping(Umin,Umax,l_min,l_max)

% perform_surftensor_remapping - compute tensor from tangent plane vectors.
%
%   T = perform_surftensor_remapping(Umin,Umax,l_min,l_max);
%
%   T is the tensor
%   Umin,Umax are eigenvector fields.
%   Cmin,Cmax are eigenvalues field.

l_min = l_min(:)';
l_max = l_max(:)';

T(1,:) = l_max.*Umax(1,:).^2 + l_min.*Umin(1,:).^2;
T(2,:) = l_max.*Umax(2,:).^2 + l_min.*Umin(2,:).^2;
T(3,:) = l_max.*Umax(3,:).^2 + l_min.*Umin(3,:).^2;
T(4,:) = l_max.*Umax(1,:).*Umax(2,:) + l_min.*Umin(1,:).*Umin(2,:);
T(5,:) = l_max.*Umax(2,:).*Umax(3,:) + l_min.*Umin(2,:).*Umin(3,:);
T(6,:) = l_max.*Umax(3,:).*Umax(1,:) + l_min.*Umin(3,:).*Umin(1,:);