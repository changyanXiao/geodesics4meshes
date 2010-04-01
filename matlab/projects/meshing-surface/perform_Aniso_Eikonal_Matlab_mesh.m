function [U,V] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_points, options)

%   perform_Aniso_Eikonal_Matlab_mesh - launch the Gauss-Seidel Iterations
%                           to compute anisotropic geodesic distances on mesh
%
%   perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, faces, T,
%   start_points, options)
%
%
%   vertex: a 3D mesh
%   start_points(i) is the index of the ith starting point.
%
%   U is the distance function to the set of starting points.
%   V is the index of the closest source point. It provide a Voronoi
%   decomposition of the domain.
%   Optional:
%   - You can provide special Initial values of the geodesic distance and
%   voronoi map, so the re-computation is done only near source points:
%       'options.U_ini' : initial U values
%       'options.V_ini' : initial V values
%
%   Copyright (c) 2010 Fethallah Benmansour


options.null = 0;

U_ini = getoptions(options, 'U_ini', []);
V_ini = getoptions(options, 'V_ini', []);

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end
start_points = start_points(:);

doUpdate = true(size(vertex));
doUpdate = getoptions(options, 'doUpdate', doUpdate);

% use fast C-coded version if possible
if exist('perform_front_propagation_2d', 'file')~=0
    if(isempty(U_ini))
        [U,V] = AnisoEikonalSolverMatlabMesh(vertex, Connectivity, T, start_points-1, doUpdate);
        V = V+1;
    else
        V_ini = V_ini-1;
        AnisoEikonalSolverMatlabMesh(vertex, Connectivity, T, start_points-1, ...
                                doUpdate, U_ini, V_ini);
        V = V_ini+1;
        U = U_ini;
    end
else
    error('You have to run compiler_mex .');
end
