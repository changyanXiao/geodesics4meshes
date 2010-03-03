% function [U, Grad] = fm2dSubGradient(h, P, w, source_point, end_points, Ob)
% Compute geodesic distance and an element of the sub-gradient with respect to the metric
% from one source point to a list of end points to reach.
%
% [U, Grad] = fm2dSubGradient(h, P, w, source_point, end_points, Ob)
% 
% P            : is the metric potential of size [Nx, Ny], should be of type double.
% h            : is the image spacing, for example h = [1/Nx;1/Ny];
% w            : is a constant real value  called regularization parameter such that P + w > 0
%                if P > 0, just take w = 0.
% source_point : is the source point in the image should be of size [2,1].
% end_points   : set of endpoint in the image domain, should be of size [2,nb_end_points]
% Ob           : (optional) add an Obstacle in the image domain. should be
%                of type logical, and of size [Nx+2, Ny+2] where the corresponding image
%                domain is (2:Nx+1, 2:Nx+1)
%
% Example      : see testGrad0* in the main directory.
%
% For theoretical and algorithmic details see:
% Derivatives with Respect to Metrics and Applications: Subgradient Marching Algorithm.
% Fethallah Benmansour, Guillaume Carlier, Gabriel Peyré and Filippo Santambrogio.
% Preprint
%
% Copyright (c) F. Benmansour 02/2010
% http://cvlab.epfl.ch/~fbenmans/