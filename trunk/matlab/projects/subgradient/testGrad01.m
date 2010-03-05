% test for fast marching
clear all; close all; clc;
%--------------------------------------------------------------------------
library_directory = '../mex/';
addpath(library_directory);
matlab_directory   = '../../matlab/';
addpath(matlab_directory);
%%
Nx = 200;
Ny = 200;
P = ones(Nx, Ny);
load 'ColorMapJet.mat';
%--------------------------------------------------------------------------
w = 1;
h = [1/Nx;1/Ny];
[source_points] = [22;22];%pick_start_points(P);
[end_points]    = [180;180];%pick_start_points(P);
%%
%--------------------------------------------------------------------------
tic
[U, Grad] = fm2dSubGradient(h, P, w, source_points, end_points);
toc
%--------------------------------------------------------------------------
figure(1); imshow(U,  [mmin(U) mmax(U)], 'InitialMagnification', 200);colormap(1-gray);
%%
figure;
imshow((Grad),  [], 'InitialMagnification', 200);colormap((1-gray).^3);
hold on
contour(Grad, 18, 'color', [0.1 0.1 0.1]);
plot(22, 22, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize', 8)
plot(180, 180, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize', 8)
axis on;
set(gca,'XTickLabel',{''})
set(gca,'XTick',[]);
set(gca,'YTickLabel',{''})
set(gca,'YTick',[]);
colorbar('vertical');
set(gca, 'fontsize', 15)
text('Interpreter','latex','String','$x_s$','Position',[20 10],'FontSize',40, 'Color', 'k');
text('Interpreter','latex','String','$x_t$','Position',[180 190],'FontSize',40, 'Color', 'k');