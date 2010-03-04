function [T] = compute_curvature_tensor_mesh(vertex, faces, Aniso)
%
%   compute_curvature_tensor_mesh - compute a curvature tensor on a
%   triangulated mesh, given a maximal anisotropy ratio
% 
%   Inputs :
%   - vertex
%   - faces
%   - Aniso : maximal anisotropy ratio
%   Output:
%   - T = exp(alpha Cmax) Umax Umax^t + exp(alpha Cmin) Umin Umin^t
%   where:
%   Umin is the direction of minimum curvature
%   Umax is the direction of maximum curvature
%   Cmin is the minimum curvature
%   Cmax is the maximum curvature
%   alpha is such that Aniso = max(sqrt(exp(alpha Cmax) / exp(alpha Cmin))) 
%   Therefore:
%   $\alpha = \frac{2 \log(Aniso)}{max(Cmax - Cmin)}
%
%   The algorithm for computing the curvatures is detailed in 
%       David Cohen-Steiner and Jean-Marie Morvan. 
%       Restricted Delaunay triangulations and normal cycle. 
%       In Proc. 19th Annual ACM Symposium on Computational Geometry, 
%       pages 237-246, 2003. 
%   and also in
%       Pierre Alliez, David Cohen-Steiner, Olivier Devillers, Bruno Le?vy, and Mathieu Desbrun. 
%       Anisotropic Polygonal Remeshing. 
%       ACM Transactions on Graphics, 2003. 
%       Note: SIGGRAPH '2003 Conference Proceedings
%
%   Copyright (c) 2010 Fethallah Benmansour --CVLab-EPFL--

T = zeros( [6, length(vertex)] );

options.verb = 0;
[Umin,Umax,Cmin,Cmax] = compute_curvature(vertex,faces,options);

alpha = 2*log(Aniso) / max(Cmax-Cmin);
l_max = exp(alpha*Cmax');
l_min = exp(alpha*Cmin');

T(1,:) = l_max.*Umax(1,:).^2 + l_min.*Umin(1,:).^2;
T(2,:) = l_max.*Umax(2,:).^2 + l_min.*Umin(2,:).^2;
T(3,:) = l_max.*Umax(3,:).^2 + l_min.*Umin(3,:).^2;
T(4,:) = l_max.*Umax(1,:).*Umax(2,:) + l_min.*Umin(1,:).*Umin(2,:);
T(5,:) = l_max.*Umax(2,:).*Umax(3,:) + l_min.*Umin(2,:).*Umin(3,:);
T(6,:) = l_max.*Umax(3,:).*Umax(1,:) + l_min.*Umin(3,:).*Umin(1,:);

return;