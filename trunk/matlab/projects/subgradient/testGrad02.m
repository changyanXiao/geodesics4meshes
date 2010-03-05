% test for fast marching
clear all; close all; clc;
%--------------------------------------------------------------------------
library_directory = '../mex/';
addpath(library_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
%%
N = 201;
[X, Y] = meshgrid(1:N, 1:N);
Sigma = 100;
P = 1./ (1.5 - exp(- sqrt((X-101).^2 + (Y-101).^2)/Sigma));
%--------------------------------------------------------------------------
w = 0;
h = [1;1]; % Or h = [1/2;1];% :-) 
Ob = zeros(N+2, N+2);
Ob = logical(Ob);
[source_points] = [22;22];%pick_start_points(P);
[end_points]    = [180;180];%pick_start_points(P);
%--------------------------------------------------------------------------
tic
[U, Grad] = fm2dSubGradient(h, P, w, source_points, end_points, Ob);
toc
%--------------------------------------------------------------------------
figure(1); imshow(U,  [mmin(U) mmax(U)], 'InitialMagnification', 200);colormap(1-gray);
%%
figure;
imshow((Grad),  [], 'InitialMagnification', 200);colormap((1-gray).^3);
hold on
contour(Grad, 25, 'color', [0.1 0.1 0.1]);
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
%%
figure; imshow(P, [], 'InitialMagnification', 200);colormap(1-gray)
colorbar('vertical');
set(gca, 'fontsize', 15)
axis on;
set(gca,'XTickLabel',{''})
set(gca,'XTick',[]);
set(gca,'YTickLabel',{''})
set(gca,'YTick',[]);

