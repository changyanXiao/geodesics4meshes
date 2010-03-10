function [U,V,Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, start_points, options)

% perform_Aniso_Eikonal_Solver_mesh - launch the Gauss-Seidel Iterations        
%                           to compute anisotropic geodesic distances on mesh
%
%   [U,V,Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T,
%   start_points, options)
%
%   Calls : number of calls of the function
%   If Calls == 0 then mesh connectivity computed
%   Else          use the last mesh connectivity
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
%       'options.U_ini_seeds' : initial U values at seed points
%       'options.V_ini_seeds' : initial V values at seed points
%   - You can provide special Initial values of the geodesic distance and
%   voronoi map, so the re-computation is done only near source points:
%       'options.U_ini' : initial U values
%       'options.V_ini' : initial V values
%
%   Copyright (c) 2010 Fethallah Benmansour & Gabriel Peyré


options.null = 0;
U_ini_seeds = getoptions(options, 'U_ini_seeds', []);
V_ini_seeds = getoptions(options, 'V_ini_seeds', []);

U_ini = getoptions(options, 'U_ini', []);
V_ini = getoptions(options, 'V_ini', []);

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end
if size(faces,1)>size(faces,2)
    faces = faces';
end
start_points = start_points(:);

doUpdate = true(size(vertex));
doUpdate = getoptions(options, 'doUpdate', doUpdate);

% use fast C-coded version if possible
if exist('perform_front_propagation_2d', 'file')~=0
    if(isempty(U_ini))
        [U,V] = AnisoEikonalSolverMesh(Calls, vertex-1, faces-1, T, start_points-1,...
                                       U_ini_seeds, V_ini_seeds, doUpdate);
        V = V+1;
        Calls = Calls + 1;
    else
        V_ini = V_ini-1;
        AnisoEikonalSolverMesh(Calls, vertex-1, faces-1, T, start_points-1, ...
                               U_ini_seeds, V_ini_seeds, doUpdate, U_ini, V_ini);
        V = V_ini+1;
        U = U_ini;
        Calls = Calls + 1;
    end
else
    error('You have to run compiler_mex before.');
end
