function path = extract_path_3d_saddle(A,saddlepoints,dUx, dUy,dUz, options)

% extract_path_2d_saddle - extract the shortest path using 
%   two gradient descents from sadlepoints.s1 and sadlepoints.s1.
%
%   path = extract_path_2d(A,end_points,options);
%
%   'A' is the distance function.
%   'end_points' is starting point (should be integer). 
%
k = length(saddlepoints);
SP = zeros(3, 2*k);
for i=1:k
  SP(:,2*i-1) = (saddlepoints(i).s1);
  SP(:,2*i) = (saddlepoints(i).s2);
end;
path = extract_path_3d(A, SP, dUx, dUy, dUz);