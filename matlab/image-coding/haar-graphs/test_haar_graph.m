% test for haar transform on graph

n = 2000;
vertex = rand(2,n);
face = compute_delaunay(vertex);
A = triangulation2adjacency(face);

v = vertex(1,:).^2 + vertex(2,:);
v = v(:);

% forward transform
vwav = perform_haar_graph(v, A, +1);

% backward transform
v1 = perform_haar_graph(vwav, A, -1);
% Check for orthogonality
err = norm(vwav)/norm(v)-1;
disp(['Should be 0: ' num2str(err) '.']);
% Check for reconstruction
err = norm(v-v1)/norm(v);
disp(['Should be 0: ' num2str(err) '.']);

% backward transform
vwavT = perform_thresholding(vwav, round(.1*n), 'largest');
[v1,w] = perform_haar_graph(vwavT, A, -1);


% original signal on an image
M0 = griddata_arbitrary(face, vertex, v, 256);
M1 = griddata_arbitrary(face, vertex, v1, 256);

clf;
imageplot({M0 M1},{'Original' 'Approximated'});

wtgt = n./[8 16 32 64];
clf;
for i=1:length(wtgt)
    % display a typical haar-let
    [tmp,k] = min(abs(w-wtgt(i)));
    vwav = zeros(n,1); vwav(k) = 1;
    v1 = perform_haar_graph(vwav, A, -1);
    Mdisp = griddata_arbitrary(face, vertex, v1, 256);
    subplot(2,2,i);
    imagesc(Mdisp);
end

