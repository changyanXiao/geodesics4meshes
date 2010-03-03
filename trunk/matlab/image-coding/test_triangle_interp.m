% test for triangulation interpolation and approximation

n = 512;
M = load_image('lena', n);
M = rescale(M);

% sample random triangles
p = 2600;

% 2600 -> 24dB alpha=1
% ->

type = 'rand';
type = 'uniform';

switch type
    case 'rand'
        vertex = floor( rand(2,p)*(n-1) + 1 );
    case 'uniform'
        n1 = floor(sqrt(p));
        n2 = round( p/n1 );
        [Y,X] = meshgrid(linspace(1,n,n2),linspace(1,n,n1));
        vertex = [X(:)'; Y(:)'];
        vertex = vertex + randn(2,p)*1e-6;
end

face = compute_delaunay(vertex);

% clf; plot_mesh(vertex,face);

% interpolation
options.interpolate = 1;
v1 = compute_orthoproj_triangulation(vertex, face, M, options);
% approximation
options.interpolate = 0;
v2 = compute_orthoproj_triangulation(vertex, face, M, options);

options.remove_nan = 1;
M1 = griddata_arbitrary(face,vertex,v1,n, options);
M2 = griddata_arbitrary(face,vertex,v2,n, options);

% error
s1 = psnr(M,M1,1);
s2 = psnr(M,M2,1);

clf;
imageplot({M M1 clamp(M2)}, {'Original', ['Interp, PSNR=' num2str(s1)], ['Approx, PSNR=' num2str(s2)]});