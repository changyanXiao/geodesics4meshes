function [start_points] = pick_start_points(M)
% pick_start_points - pick start points for front propagation.
%
%   [start_points] = pick_start_points(M);
%   
%   The user right-click on a set of point.
%   Another left click stop the process.

clf;
hold on;
if ~isempty(M)
imshow(M, [mmin(M) mmax(M)], 'InitialMagnification', 100);
axis image;
axis off;
colormap gray(256);
end

m = size(M);

b = 1;
start_points = [];
while b(end)==1
    [y1,x1,b] = ginput(1);
    if b==1
        start_points = [start_points [x1;y1]];
    end;
end;