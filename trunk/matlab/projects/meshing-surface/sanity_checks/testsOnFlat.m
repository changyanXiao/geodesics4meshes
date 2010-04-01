%% Sanity Check of anisotropic geodesic distances with the new struture
clear all; close all; clc
add_base_paths;
%%
load 'FlatMesh';
nverts = length(vertex);
%% chose the anisotropy
T = zeros(6, nverts);
theta = -pi/6;
alpha = 3;
beta = 1;
e_theta = [cos(theta), sin(theta)];
e_theta_T = [-sin(theta), cos(theta)];
tensor = alpha^2 *e_theta'*e_theta + beta^2*e_theta_T'*e_theta_T;
T(1,:) = tensor(1, 1);
T(2,:) = tensor(2, 2);
T(4,:) = tensor(1, 2);
%% start from both sources
start_points =[7432,3923];%, 100, 10000, 1];
options.start_points = start_points;
tic
[U, V] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_points, options);
toc
%% compute real distance and real Voronoi
U_th = 1e9*ones(size(U));
V_th = zeros(size(U));
for k = 1:length(start_points)
    dist = sqrt(T(1,:).*(vertex(1, :) - vertex(1, start_points(k))).^2 +...
                T(2,:).*(vertex(2, :) - vertex(2, start_points(k))).^2 +...
                T(3,:).*(vertex(3, :) - vertex(3, start_points(k))).^2 +...
                2*T(4,:).*(vertex(1, :) - vertex(1, start_points(k))).*(vertex(2, :) - vertex(2, start_points(k)))+...
                2*T(5,:).*(vertex(2, :) - vertex(2, start_points(k))).*(vertex(3, :) - vertex(3, start_points(k)))+...
                2*T(6,:).*(vertex(3, :) - vertex(3, start_points(k))).*(vertex(1, :) - vertex(1, start_points(k))));
	V_th(dist' <= U_th) = k;
	U_th = min(U_th , dist');
end
%%
error = abs(U_th -U);
U_th(start_points) = 1;
relative_error = error./U_th;
fprintf('start from both sources, max relative error is %e\n', max(relative_error));
%%
figure;
plot_fast_marching_mesh(vertex,faces, relative_error, [], options);
zoom(.7);
title('start from both sources: relative distance error');
%%
figure;
plot_fast_marching_mesh(vertex,faces, V_th, [], options);
zoom(.7);
title('theoretical voronoi');
%%
figure;
plot_fast_marching_mesh(vertex,faces, double(V)-V_th, [], options);
zoom(.7);
title('start from both sources: voronoi error');
%% start from each source and take the minimum
start_point1 =7432;
tic
[U1, V1] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_point1, options);
toc
start_point2 =3923;
tic
[U2, V2] = perform_Aniso_Eikonal_Matlab_mesh(vertex, Connectivity, T, start_point2, options);
toc
UU = min(U1, U2);
V(UU==U1) = 1;
V(UU==U2) = 2;
%% compute real distance and real Voronoi
U_th = 1e9*ones(size(U));
V_th = zeros(size(U));
for k = 1:length(start_points)
    dist = sqrt(T(1,:).*(vertex(1, :) - vertex(1, start_points(k))).^2 +...
                T(2,:).*(vertex(2, :) - vertex(2, start_points(k))).^2 +...
                T(3,:).*(vertex(3, :) - vertex(3, start_points(k))).^2 +...
                2*T(4,:).*(vertex(1, :) - vertex(1, start_points(k))).*(vertex(2, :) - vertex(2, start_points(k)))+...
                2*T(5,:).*(vertex(2, :) - vertex(2, start_points(k))).*(vertex(3, :) - vertex(3, start_points(k)))+...
                2*T(6,:).*(vertex(3, :) - vertex(3, start_points(k))).*(vertex(1, :) - vertex(1, start_points(k))));
	V_th(dist' < U_th) = k;
	U_th = min(U_th , dist');
end
%%
error = abs(U_th -UU);
U_th(start_points) = 1;
relative_error = error./U_th;
fprintf('start from each source and take the minimum, max relative error is %e\n', max(relative_error));
%%
figure;
plot_fast_marching_mesh(vertex,faces, relative_error, [], options);
zoom(.7);
title('start from each source and take the minimum: relative distance error');
%%
figure;
plot_fast_marching_mesh(vertex,faces, double(V)-V_th, [], options);
zoom(.7);
title('start from each source and take the minimum:voronoi error');
%%
figure;
plot_fast_marching_mesh(vertex,faces, abs(UU-U), [], options);
zoom(.7);
title('Diff computed distances');