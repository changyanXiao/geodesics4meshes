% test for error decay with M

n = 256;
M = load_image('lena', n);
M = rescale(M);

% sample random triangles
plist = round(linspace(300,5000,10));
ntests = length(plist);

err_interp = [];
err_approx = [];
vertex = [];

for i=1:ntests
    p = plist(i);
    dp = p - size(vertex,2);
    vertex(:,end+1:p) = floor( rand(2,dp)*(n-1) + 1 );

    face = delaunay(vertex(1,:), vertex(2,:))';


    % interpolation
    v1 = interp2(M, vertex(2,:), vertex(1,:));
    % approximation
    v2 = compute_orthoproj_triangulation(vertex, face, M);

    options.remove_nan = 1;
    M1 = griddata_arbitrary(face,vertex,v1,n, options);
    M2 = griddata_arbitrary(face,vertex,v2,n, options);

    % error
    err_interp(end+1) = snr(M,M1);
    err_approx(end+1) = snr(M,M2);
end

clf;
plot(plist,err_interp, plist,err_approx);
xlabel('Nbr.Vertex'); ylabel('SNR');
legend('Interpolation', 'Approximation');