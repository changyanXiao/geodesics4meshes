clear all; close all; clc;
%%
library_directory = '../mex';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
meshes_directory   = '../meshes/';
addpath(meshes_directory);
%%
Nb_calls = 0;
clear options;
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 0;
[vertex,faces] = compute_semiregular_sphere(7,options);
nverts = size(vertex,2);
nstart = 1;
start_points = [];
%%
Calls = 0;
%% D'abord sans anisotropie
T = zeros(6, nverts);
T(1,:) = 1;
T(2,:) = 1;
T(3,:) = 1;
%%
figure;
added_points = 1;
NB_points = 6;
start_points = [ start_points floor(rand(nstart,1)*nverts)+1];
options.start_points = start_points;
tic
[U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, start_points);
toc
while added_points < NB_points
    options.start_points = start_points;
    plot_fast_marching_mesh(vertex,faces, V, [], options);
    zoom(.7);
    pause;
    options.U_ini = U;
    options.V_ini = V;
    pointsMax = find(U == max(U(:)));% langth(pointsMax) >= 1, 
    start_points = [start_points pointsMax(1)];% therfore, take only one
    tic
    [U, V, Calls] = perform_Aniso_Eikonal_Solver_mesh(Calls, vertex, faces, T, start_points(end), options);
    toc
    added_points = added_points+1;
end
%%
options.start_points = start_points;
figure;
plot_fast_marching_mesh(vertex,faces, U, [], options);
zoom(.7);
%%
figure;
plot_fast_marching_mesh(vertex,faces, V, [], options);
zoom(.7);