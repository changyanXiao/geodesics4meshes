function [UU, Calls] = perform_geodesic_computation_extended(Calls, vertex, faces, W, T, U, V, Voronoi_index, Vext)

% perform_geodesic_computation_extended - 
%
% [UU, calls] = perform_geodesic_computation_extended(Calls, vertex, faces,
% W, T, U, V, Voronoi_index, Vext)
%
%   Calls : number of calls of the function
%   If Calls == 0 then mesh connectivity computed
%   Else          use the last mesh connectivity
%
%   vertex, faces, W: a 3D mesh, W being the adjacency matrix
%   T : the metric tensor
%   U and V are the initial geodesic distance and associated voronoi.
%   Voronoi_index: is the index of the voronoi regio to extend
%   Vext         : the width of the extension
%   
%   returns:
%   UU : geodesic distance on the extended domain
%   Calls is incremented
%
%   Copyright (c) 2010 Fethallah Benmansour & Gabriel PeyrŽ



v = double(V==Voronoi_index);
for iext=1:Vext
	v = W*v;
end
v = v>0;
% perform propagation on the enlarged region
options.doUpdate = v;
U_ini = U;
U_ini(V~=Voronoi_index) = 1e6;% give the extended part a large value (infinity)
V_ini = int16(Voronoi_index*ones(size(vertex)));
V_ini(V~=Voronoi_index) = Voronoi_index; % give the extension the same voronoi
options.U_ini = U_ini;
options.V_ini = V_ini;
% pour forcer l'undate seulemnt dans l'extension, il suffit de donner un
% seed sur le bord
seed = find( (U == max(U(V==Voronoi_index))) & (V == Voronoi_index) );
U_ini_seeds = U(seed);
V_ini_seeds = V(seed);
options.U_ini_seeds = U_ini_seeds;
options.V_ini_seeds = V_ini_seeds;

tic
[UU, V1, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, seed, options);
toc
UU(~v) = 0;
clear V1;
return;