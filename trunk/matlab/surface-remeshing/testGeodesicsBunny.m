clear all; close all; clc;
%%
library_directory = '../mex';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
meshes_directory   = '../meshes/';
addpath(meshes_directory);
%%
[vertex,faces] = read_mesh('bunny.off');
nvert = size(vertex,2);nverts = size(vertex,2);
nstart = 1;
%%
%start_points = floor(rand(nstart,1)*nverts)+1;
start_points = 160;
options.start_points = start_points;
%% 
Aniso = 10000;
[T] = compute_curvature_tensor_mesh(vertex, faces, Aniso);
%%
tic
[U, V, dUx, dUy, dUz] = perform_Aniso_Eikonal_Solver_mesh(vertex, faces, T, start_points);
toc
%%
figure;
plot_fast_marching_mesh(vertex,faces, U2, [], options);
zoom(.7);
%%
figure;
plot_fast_marching_mesh(vertex,faces, V, [], options);
zoom(.7);