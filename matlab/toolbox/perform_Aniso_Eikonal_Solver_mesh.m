function [U,V,dUx,dUy,dUz] = perform_Aniso_Eikonal_Solver_mesh(vertex, faces, T, start_points, options)

% perform_Aniso_Eikonal_Solver_mesh - launch the Gauss-Seidel Iterations        
%                           to compute anisotropic geodesic distances on mesh
%
%   [U,V,dUx,dUy,dUz] = perform_Aniso_Eikonal_Solver_mesh(vertex, faces, T,
%   start_points, options)
%
%   vertex, faces: a 3D mesh
%   start_points(i) is the index of the ith starting point.
%
%   U is the distance function to the set of starting points.
%   V is the index of the closest source point. It provide a Voronoi
%   decomposition of the domain.
%   dUx, dUy and dUz provides the direction of characteristics (namely geodesics directions)
%   Optional:
%   - You can provide special Initial values of the geodesic distance at
%   the source points conditions for stop in options : 
%       'options.U_ini' : initial U values
%       'options.V_ini' : initial V values
%
%   Copyright (c) 2010 Fethallah Benmansour & Gabriel Peyré


options.null = 0;
U_ini = getoptions(options, 'U_ini', []);
V_ini = getoptions(options, 'V_ini', []);

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end
if size(faces,1)>size(faces,2)
    faces = faces';
end
start_points = start_points(:);

% use fast C-coded version if possible
if exist('perform_front_propagation_2d', 'file')~=0
    [U,V,dUx,dUy,dUz] = AnisoEikonalSolverMesh(vertex-1, faces-1, T,start_points-1,U_ini, V_ini);
    V = V+1;
else
    error('You have to run compiler_mex before.');
end
