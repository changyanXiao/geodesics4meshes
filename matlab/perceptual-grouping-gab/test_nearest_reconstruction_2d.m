%  test for 2D contour reconstruction with potential
%
%  Copyright (c) 2008 Sébastien Bougleux, Gabriel Peyré

if 1
    path(path, 'data/');
    name = 'curves';
    n = 150;
    M = load_image(name, n);
    M = rescale(clamp(M,0,255));
    W = 1 - M + 0.00001;

    rep = 'results/recons-2d/';
    if not(exist(rep))
        mkdir(rep);
    end

    warning off;
    imwrite(rescale(W), [rep name '-map.png'], 'png');
    warning on;
end
%W = 1 - M + 0.00001;
options.nb_iter_max = Inf;
disp('Performing front propagation.');
load points1;
[D,S,Q,A,NN] = perform_contour_completion(W, start_points, options);
npaths=size(NN,2);

% plot distance, saddle points and topological graph
figure;
D1 = perform_histogram_equalization(D, linspace(0,1,n^2));
imagesc(D1);
axis image; axis off;
colormap jet(256);
hold on;
for i=1:nstart
    for j=1:nstart
        if ((i>j) & (A(i,j) ~= 0.))
            h = line([start_points(2,i),start_points(2,j)],[start_points(1,i),start_points(1,j)], 'color', 'k');
            set(h,'LineWidth',3);
        end
    end
end
h = plot(start_points(2,:),start_points(1,:), 'k.');
set(h, 'MarkerSize', 30);
j = 1;
stop = npaths/2;
for i=1:stop
    h = plot((NN(2,j)+NN(2,j+1))/2.,(NN(1,j)+NN(1,j+1))/2., 'k+');
    set(h, 'MarkerSize', 12);
    set(h, 'LineWidth', 4);
    j=j+2;
end
saveas(gcf, [rep name '-dist.eps'], 'psc2');
hold off;

% paths
figure;
imagesc(M);
axis image; axis off;
colormap gray(256);
hold on;
for i=1:npaths
    %A=D;
    %B=find(Q~=Q(end_points(1,i),end_points(2,i)));
    %A(B)=Inf;
    paths{i} = compute_geodesic(D,NN(:,i));
    h = plot( paths{i}(2,:), paths{i}(1,:), 'r' );
    set(h, 'LineWidth', 3);    
end
j=1;
for i=1:stop
    h = plot((NN(2,j)+NN(2,j+1))/2.,(NN(1,j)+NN(1,j+1))/2., 'b+');
    set(h, 'MarkerSize', 12);
    set(h, 'LineWidth', 4);
    j=j+2;
end
h = plot(start_points(2,:),start_points(1,:), 'b.');    
set(h, 'MarkerSize', 30);
saveas(gcf, [rep name '-recons.eps'], 'psc2');
hold off;
