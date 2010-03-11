function [flist,lambda] = compute_triple_points(vertex,faces, landmarks, D)

% compute_triple_points - compute the location of triple points
%
%   [flist,lambda] = compute_triple_points(vertex,face, landmarks, ULoc);
%
%   ULoc(:,i) is the distance to landmarks(i)
%
%   flist is the list of face containing triple points
%   l=lambda(:,i) is the barycentric coordinate in f=flist(i,:)
%   so that the position of the triple point is 
%       l .* vertex(:,f)
%
%   Copyright (c) 2010 Gabriel Peyre

n = size(vertex,2);

%%
% Compute the 3 closest voronoi indexes at each point.

D1 = D;
Q = []; DQ = [];
for i=1:3
    [DQ(:,i),Q(:,i)] = min(D1,[],2);
    Qind = (1:n)' + (Q(:,i)-1)*n;
    D1(Qind) = Inf;
end
Q1 = Q(:,1); % first voronoi

%% 
% Find faces with triple points.

A = (Q1(faces(1,:))~=Q1(faces(2,:))) + (Q1(faces(2,:))~=Q1(faces(3,:))) + (Q1(faces(1,:))~=Q1(faces(3,:)));
I = find(A==3);
flist = faces(:,I);
lambda = [];
for i=I(:)'
    % face
    f = faces(:,i);  
    % voronoi indexes
    v = Q1(f);
    % corresponding distances
    d = D( f,v );
    % d(:,k) is the distance for the kth function,
    B = [ d(:,1)'-d(:,2)'; d(:,1)'-d(:,3)'; [1 1 1]];
    % B(k,:) gives a distance at points x
    lambda(:,end+1) = B \ [0;0;1];
end
% enforce within triangle property
if not(isempty(lambda))
    lambda = max(lambda,0);
    lambda = lambda ./ repmat(sum(lambda), [3 1]);
end