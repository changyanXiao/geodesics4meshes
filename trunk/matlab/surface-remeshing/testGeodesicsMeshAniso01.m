clear all; close all; clc;
%%
library_directory = '../mex';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
meshes_directory   = '../meshes/';
addpath(meshes_directory);
%%
clear options;
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 0;
[vertex,faces] = compute_semiregular_sphere(7,options);
nverts = size(vertex,2);
nstart = 10;
%%
start_points = floor(rand(nstart,1)*nverts)+1;
options.start_points = start_points;
%% D'abord sans anisotropie
T = zeros(6, nverts);
T(1,:) = 1;
T(2,:) = 1;
T(3,:) = 1;
%%
tic
[U, V, dUx, dUy, dUz] = perform_Aniso_Eikonal_Solver_mesh(vertex, faces, T, start_points);
toc
%%
figure;
plot_fast_marching_mesh(vertex,faces, U, [], options);
zoom(.7);
%%
figure;
plot_fast_marching_mesh(vertex,faces, V, [], options);
zoom(.7);