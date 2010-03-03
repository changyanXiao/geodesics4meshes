% test for metric uniformization
n = 150;

% obstacles
x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);
r = 0.2;
c = [0.5 0.5];
R = (X-c(1)).^2 + (Y-c(2)).^2 < r^2;
% constraint map
L = zeros(n) + Inf;
L(R==1) = -Inf;
options.constraint_map = L;

% start / end
clf;
imagesc(R); axis image;
[y,x] = ginput(2);
x0 = round([x(1);y(1)]);
x1 = round([x(2);y(2)]);

% initial potential
W = ones(n);
% W(R==0) = W(R==0)/sum(W(R==0));
W(R==1) = 0;

% iterate
K = 8*n; % number of points updated each time
niter = 500;
for i=1:niter
    progressbar(i,niter);
    % perform propagation
    [D0,S] = perform_fast_marching(W, x0, options);
    [D1,S] = perform_fast_marching(W, x1, options);
    % distance between the two points
%    d = ( D0(x1(1),x1(2)) + D1(x0(1),x0(2)) )/2;   
    d = min( D0(x1(1),x1(2)), D1(x0(1),x0(2)) );
    % map to threshold
    w = (D0+D1-d);
    % find threshold
    s = sort(w(R==0));
    epsilon = s(K);
    I = find(w<=epsilon);
    imagesc(w<=epsilon);
    % update metric
    % W(I) = W(I)/1.5;
    rho = 1/10;
    W = W*(1-rho/4);
    W(I) = W(I) + rho;
    % rescale
%    W(R==0) = W(R==0)/sum(W(R==0));
    W(R==1) = 0;
end
U = 1./(0.00001+W);
U(R==1) = 0;
imagesc(U);
