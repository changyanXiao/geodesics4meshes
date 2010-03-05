function [D,Q] = perform_distance_transform(mask, points)

% [D,Q] = perform_distance_transform(mask, points);

I = find(mask);
[x,y] = ind2sub(size(mask),I);
base = [x(:) y(:)]';

D1 = compute_distance_to_points(points,base);
[D1,Q1] = min(D1, [], 2);
D = mask*0;
Q = mask*0;
D(I) = sqrt(D1);
Q(I) = Q1;

return;

D = mask; Q = mask;
D(mask==0) = Inf;
Q(mask==0) = Inf;
D(I) = reshape( sqrt(D1), size(mask));
Q(I) = reshape(Q1, size(mask));
