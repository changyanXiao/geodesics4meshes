%% Sanity Check of anisotropic geodesic distances with the new struture
clear all; close all; clc
add_base_paths;
%%
Trig = LoadTrig('flat');
vertex = Trig.Coords';
nverts = length(vertex);
faces = Trig.Facets';
%% compute connectivity
tic
Connectivity = ComputeMeshConnectivity(vertex-1, faces-1);
toc
%% D'abord sans anisotropie
T = zeros(6, nverts);
T(1,:) = 1;
T(2,:) = 1;
T(3,:) = 1;
%%
start_points =[7432,3923];
options.start_points = start_points;
tic
[U, V] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_points, options);
toc
%%
% start_points =[7432];
% %start_points = [1,10000];
% options.start_points = start_points;
% tic
% [U1, V1] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_points, options);
% toc
% %%
% start_points =[3923];
% %start_points = [1,10000];
% options.start_points = start_points;
% tic
% [U2, V2] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_points, options);
% toc
%%
figure;
plot_fast_marching_mesh(vertex,faces, min(U1,U2)-U, [], options);
zoom(.7);
%%
figure;
plot_fast_marching_mesh(vertex,faces, V, [], options);
zoom(.7);