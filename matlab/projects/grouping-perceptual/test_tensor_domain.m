% compute a tensor field inside a shape

rep = ['results/tensor-domain/'];
if not(exist(rep))
    mkdir(rep);
end


n = 200;
b = 1;
bound = [];
while b==1
    clf;
    if not(isempty(bound))
    h = plot(bound(2,:), bound(1,:), '.-');
    set(h, 'LineWidth', 2);
    end
    title('Click [left] for new point, [right] to exit');
    axis([1 n 1 n]); axis on; box on;
    [y,x,b] = ginput(1);
    bound(:,end+1) = [x;y];
end

bound(:,end+1) = bound(:,1);

% digitalize the boundary shape
sk = draw_polygons(n,1,bound);

% compute the inside
L = ones(n)+Inf; L(sk==1) = -Inf;
options.constraint_map = L;
[D,S,Q] = perform_fast_marching(ones(n), [1 1]', options);
shape = (S==1);

clf;
hold on;
imageplot(shape); axis image;
h = plot(bound(2,:), bound(1,:), '.-'); set(h, 'LineWidth', 2);


% compute a tensor field
b = 1;
mask = zeros(n);
T = zeros(n,n,3);
while b==1
    for k=1:2
        title(['Click [left] for point #' num2str(k) ', [right] to exit']);clf;
        hold on;
        imageplot(shape+mask); axis image;
        h = plot(bound(2,:), bound(1,:), '.-'); set(h, 'LineWidth', 2);
        [y(k),x(k),b] = ginput(1);
        if b~=1
            break
        end
    end
    % compute the tensor field
    v = [diff(x);diff(y)]; v = v/norm(v, 'fro');
    T0 = v*v';
    
    % position of points along the line
    p = 200;
    t = linspace(0,1, p);
    A = repmat(t,[2,1]) .* repmat( [x(1);y(1)],[1 p]) + ... 
        repmat(1-t,[2,1]) .* repmat( [x(2);y(2)],[1 p]);
    A = round(A);
    I = A(1,:) + (A(2,:)-1)*n;
    mask(I) = 1;
    
    T(I) = T0(1,1);
    T(I+n^2) = T0(2,2);
    T(I+2*n^2) = T0(2,1);
        
end

%% perform diffusion
options.mask = shape;
options.niter = 600;
T1 = perform_tensor_densification(T, mask, options);

U = perform_tensor_mapping(T1,1);
U(:,:,2) = .6;
T1 = perform_tensor_mapping(U,-1);

clf;
options.sub = 4;
plot_tensor_field(T1, shape-.5*mask, options);
saveas(gcf, [rep 'tensor-domain.png']);