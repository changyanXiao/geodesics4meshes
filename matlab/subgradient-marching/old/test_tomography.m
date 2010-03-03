% test for travel time tomography

% the energy to optimize is
%   E(xi) = \sum_{s,t} | d_{xi}(s,t) - d_{s,t} |^2
% where d_xi is the geodesic distance for metric xi
% and d_{s,t} are the available distances
%
% The gradient of the energy is 
%   grad_xi(E) = \sum_{s,t} grad_xi(s,t) * ( d_{xi}(s,t) - d_{s,t} )  
% where grad_xi(s,t) is the vector gradient of the distance between s and t
% with respect to the metric xi.

path(path, 'data/');
path(path, 'toolbox/');

n = 50;
name = 'mountains';
name = 'bump';
name = 'constant';
options.center = [.4 .4];
options.bump_size = .15;
M = load_image(name, n, options);
if strcmp(name, 'mountains')
    M = -M;
end

w = 1e-3;   % regularization
h = [1;1];  % spacing
% potential one wish to find
P0 = rescale(M,.1,1);


% points on the boundary
delta = 1;
s = delta:n-delta+1;
boundary = [[s;s*0+delta] [s;s*0+n-delta+1] [s*0+delta;s] [s*0+n-delta+1;s]];
nbound = size(boundary,2);

% starting points
nstart = 60;
istart = floor( rand(nstart,1)*nbound )+1;
start_points = boundary(:,istart);
start_points = round([n;n]/2);
nstart = size(start_points,2);

% ending points
nend = nbound;
sel = randperm(nbound); sel = sort(sel(1:nend));
end_points = boundary(:,sel);


nend = n^2;
sel = randperm(n^2); sel = sort(sel(1:nend));
[x y] = ind2sub([n n], sel);
end_points = [x(:)'; y(:)'];

% indexes of ending points
Iend = end_points(1,:) + (end_points(2,:)-1)*n;

% initialize the set of distances
D = zeros(nstart,nend);
options.end_points = [];
disp('--> Collecting initial distances.');
for i=1:nstart
%    U0 = perform_fast_marching(1./P0, start_points(:,i), options)*n;
    [U0, Grad] = fm2dGradient(h, double(P0), w, start_points(:,i), end_points);
    D(i,:) = U0( Iend )';
end

% initialize the map
P = zeros(n)+mean(P0(:))/4; % P = P0;
niter = 200;
tau = .05;
Pswg = zeros(n,n,niter);
E = zeros(niter-1,1); % energy 
for it=2:niter
    progressbar(it,niter);
    % compute the gradient
    G = zeros(n);
    for i=1:nstart % iterate on starting points
        % remove point to close to start
        d = sqrt( sum( ( end_points-repmat(start_points(:,i), [1 nend]) ).^2, 1) );
        sel = find( d>=.01*n );
        % perform gradien computation
        [U, Grad] = fm2dGradient(h, double(P), w, start_points(:,i), end_points(:,sel));
        d = ( U(Iend(sel)) - D(i,sel) );
        Grad = Grad.*repmat( reshape(d,1,1,length(sel)), [n n 1] );
        % accumulate gradient
        G = G + mean(Grad,3);
        % compute the energy
        E(it-1) = E(it-1) + sum(d.^2);
    end
    % perfomr descent
    P = P - tau*G;
    % enforce smoothness
%    P = perform_blurring(P,2);
    % crop values 
    % P = clamp(P, min(P0(:)), max(P0(:)) );
    % display
    imageplot(P); drawnow;
    % to make a movie
    Pswg(:,:,it) = rescale(P);
end

   