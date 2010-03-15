function [T,Umin,Umax] = compute_surface_metric(vertex,faces, metric_type, options)

% compute_surface_metric - compute a tensor metric for meshes.
%
%   [T,Umin,Umax] = compute_surface_metric(vertex,faces, metric_type, options);
%
%   It is a curvature-based metric.
%
%   metric_type is either 'isotropic', 'anisotropic', 'constant'.
%
%   Copyright (c) 2010 Gabriel Peyre

options.verb = 0;
[Umin,Umax,Cmin,Cmax] = compute_curvature(vertex,faces,options);
Cmin = abs(Cmin); Cmax = abs(Cmax);
% ensure Cmin<Cmax
I = find(Cmin>Cmax);
Cmax1=Cmax; Cmax(I) = Cmin(I); Cmin(I) = Cmax1(I);
Umax1=Umax; Umax(:,I) = Umin(:,I); Umin(:,I) = Umax1(:,I);
% ensure strong convexity
%epsilon = .000;
%Cmin = Cmin+epsilon; Cmax = Cmax + epsilon;

A = Cmax./Cmin;
E = Cmax+Cmin;

switch metric_type
    case 'isotropic'
        A = A*0+1;
%        E = perform_hist_eq(A, 'linear');
%        E = rescale(E,1,10);
    case 'isotropic-boost'
        A = A*0+1;
        E = E.^1.5;
    case 'constant'
        A = A*0+1;
        E = E*0+1; 
    case 'anisotropic'
    case 'anisotropic-boost'
        A = A*2;
    otherwise
        error('Unknown metric type.');
end

Cmin = E./(1+A);
Cmax = Cmin.*A;

T = perform_surftensor_remapping(Umin,Umax,Cmin,Cmax);
