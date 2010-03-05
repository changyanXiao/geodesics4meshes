% test for maximization of distance
% under L1 and Linf constraints on the metric

n = 160; % size of the problem

path(path,'data/');

name = 'square';
name = 'room';
name = 'cavern';
name = 'simple';

rep = 'results/maxim-dist/';
if not(exist(rep))
    mkdir(rep);
end

% bounds on the metric, should be around 1
rmax = 1.1;
rmin = .9;
tau = .3;

% creat obstacle
switch name
    case 'simple'
        % no obstacle
        Ob = zeros(n+2);
        eta = .3;
        pstart = round(eta*[1 1]'*n); pend = round((1-eta)*[1 1]'*n);
        rmin = 0; rmax = 5;
        tau = .3;
    case 'square'
        r = .4; % size of the square
        Ob = zeros(n+2);
        x = linspace(-1,1,n+2);
        [Y,X] = meshgrid(x,x);
        Ob( max(abs(X),abs(Y))<r ) = 1;
        pstart = round([.1 .1]'*n); pend = round([.9 .9]'*n);
        rmin = .8; rmax = 1.2;
    case {'cavern' 'room'}
        Ob = rescale(load_image(name, n+2));
        Ob = Ob<.5;
        if strcmp(name, 'cavern')
            Ob = perform_blurring(Ob, 8);
            Ob = Ob<.5;
            Ob = perform_blurring(Ob, 8);
            Ob = Ob<.5;
        end
        if Ob(1)==0
            Ob = 1-Ob;
        end
        if not(exist('pend'))
        [pstart,pend] = pick_start_end_point(Ob);
        end
        rmin = .7; rmax = 1.3;
        rmin = .8; rmax = 1.5; % for room
        tau = .2;
end
Ob = logical(Ob);
Mask = (1-Ob(2:end-1,2:end-1));

% parameters for Feth code
w = 1e-3;   % regularization
h = [1;1];  % spacing

% number of available grid cells
v = sum(Mask(:));
% initial flat metric
P = ones(n).*Mask;


niter = 500;
taulist = linspace(.8, .3, niter);
Psvg = zeros(n,n,niter);
for i=1:niter   
    progressbar(i,niter);
    if mod(i,2)==1
        [U, Grad] = fm2dGradient(h, max(double(P),1e-10), w, pstart, pend,Ob);
    else
        [U, Grad] = fm2dGradient(h, max(double(P),1e-10), w, pend, pstart,Ob);
    end
    P = P + taulist(i)*Grad;
    P(Mask==1) = compute_projection_linf(P(Mask==1),rmin,rmax);
    if 0
    for k=1:5
        P(Mask==1) = max( min(P(Mask==1),rmax), rmin);
        P(Mask==1) = P(Mask==1)-mean(P(Mask==1))+1;
    end
    end
    % for the display: put back
    A = clamp(P,rmin,rmax); A(1) = rmin; A(2) = rmax;
    A = rescale(A); 
    A(Mask==0) = -0.15; A = rescale(A);    
    imageplot(A); drawnow;
    Psvg(:,:,i) = A;
end

str = ['-min=' num2str(rmin) '-max=' num2str(rmax)];
warning off;
imwrite(rescale(A), [rep name str '.png'], 'png');
warning on;
options.format = 'gif';
compute_movie_file(Psvg(:,:,1:i-1), [rep name str '.gif'], options);