
%% testing the error decay
ntests = 12;
alphalist = [1/2 1];
mlist = linspace(100, 1500, ntests);
err_geod = [];
err_eucl = [];
err_isot = [];
for i=1:ntests
    m = mlist(i);
    % approximate with anisotropic geodesic delaunay
    for k=1:length(alphalist)
        alpha = alphalist(k);
        T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);
        [U, dUx,dUy,V, vertex,faces] = FPSCAniso2D([1;1], T, m, 0); faces = faces+1;
        v1 = compute_orthoproj_triangulation(vertex, faces, M);
        M1 = griddata_arbitrary(faces,vertex,v1,n, options);
        err_geod(i,k) = snr(M, M1);
    end
    % approximate with euclidean delaunay
    faceseuc = compute_delaunay(vertex);
    v1 = compute_orthoproj_triangulation(vertex, faceseuc, M);
    M1 = griddata_arbitrary(faceseuc,vertex,v1,n, options);
    err_eucl(i) = snr(M, M1);
    % save_image(clamp({M M1}), {'original' ['m-' num2str(m)]}, options);
end

clf; hold on;
h = plot(mlist, err_geod(:,1), 'k-'); 
set(h, 'LineWidth', lw); set(h, 'MarkerSize', 1);
h = plot(mlist, err_geod(:,2), 'k.-'); 
set(h, 'LineWidth', lw); set(h, 'MarkerSize', 15);
h = plot(mlist, err_eucl(:), 'k*-'); 
set(h, 'LineWidth', lw); set(h, 'MarkerSize', 10);
axis tight;
legend( ['Anisotropic, \alpha=' num2str(alphalist(1))], ...
        ['Anisotropic, \alpha=' num2str(alphalist(2))], 'Euclidean');
saveas(gcf, [rep name '-influ-m.eps'], 'eps');