%% Test for image coding with anisotropic triangulations

%% Load the image
name = 'cameraman';
name = 'triangle';
name = 'peppers';

rep = 'results/anisotropic-coding/';
if not(exist(rep))
    mkdir(rep);
end

%% Parameter for the tensor field remapping
n = 256; n0 = n;
switch name
    case 'triangle'
        n = 512; n0 = n;
        epsilon = 1e-6;
        alpha = .7;
        m = 200;
        sigma1 = 2;   % for gradient computation
        sigma2 = 6;     % for filtering the tensor field
    case 'cameraman'
        n = 256; n0 = n;
        % for cameraman
        epsilon = 1e-6;
        alpha = .75;
        m = 500;
        sigma1 = 2;   % for gradient computation
        sigma2 = 4;     % for filtering the tensor field
    case 'peppers'
        % for peppers
        epsilon = 1e-6;
        alpha = .75;
        m = 600;
        sigma1 = 2;   % for gradient computation
        sigma2 = 6;     % for filtering the tensor field
        n = 256; n0 = 512;
end


M = load_image(name, n0);
M = rescale(crop(M,n));
% M = perform_blurring(M,3);


%% Parameter for interpolation
options.remove_nan = 1;
%% Structure tensor extraction
options.use_anisotropic = 0;
T0 = compute_structure_tensor(M, sigma1, sigma2);
[e1,e2,l1,l2] = perform_tensor_decomp(T0);
T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);

if 0
%% use hessian
H = compute_hessian( perform_blurring(M,3), options);
[e1,e2,l1,l2] = perform_tensor_decomp(H, 'abs');
T = perform_tensor_recomp(e1, e2, (abs(l1)+epsilon).^alpha, (abs(l2)+epsilon).^alpha);
T = perform_blurring(T,3);
end

% Farthest seeding
fprintf('--> Farthest point seeding ... ');
[U, dUx,dUy,V, vertex,faces] = FPSCAniso2D([1;1], T, m, 0);
m1 = size(vertex,2); % real number of vertices
fprintf('done.\n');
faces = faces+1;

% quantization of positions
Tpos = 2;
[vertexT, vertexI] = perform_quantization(vertex, Tpos, 1);

% approximated image
options.verb = 0;
v0 = compute_orthoproj_triangulation(vertexT, faces, M);
v = clamp(v0,-.2,1.2);
% quantized values
Tval = 1/20;
[vT, vI] = perform_quantization(v, Tval, 1);

if 0
% ortho change of basis
TvalU = .5;
A = triangulation2adjacency(faces);
[U,S] = eig(full(A));
vU = U'*v;
[vUT, vI] = perform_quantization(vU, TvalU, 1);
vT = U*vUT;
end

M1T = griddata_arbitrary(faces,vertexT,vT, n);

% for reference onlY ...
if 0
v1 = compute_orthoproj_triangulation(vertex, faces, M);
M1  = griddata_arbitrary(faces,vertex,v1, n);
end

%% Count bits using fixed
RposF = length(vertexI(:))*log2( mmax(vertexI) - mmin(vertexI) + 1 );
RvalF = length(vI(:)) * log2( mmax(vI)-mmin(vI)+1 );
%% count bits using entropy
RposE = length(vertexI(:))*compute_entropy(vertexI(:));
RvalE = length(vI(:))*compute_entropy(vI(:));
%% count bits using arithmetic coder
[tmp,RposA] = perform_arithmetic_coding(vertexI(:), +1, options);
[tmp,RvalA] = perform_arithmetic_coding(vI(:), +1, options);

%% Connectivity coding
nbits0 = round( (RposA+RvalA)*.8 );
Jmin = 4;
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;
MW = perform_wavelet_transform(M*255, Jmin, +1, options);
[MWT,nbits] = perform_jp2k_degradation(MW,Jmin,nbits0,M, 8);
Mwav = perform_wavelet_transform(MWT, Jmin, -1, options);
Mwav = Mwav/255;

egeod = snr(M,M1T);
ewav = snr(M,Mwav);
clf;
imageplot( clamp(M1T), ['Geodesic, SNR=' num2str(egeod,3) 'dB'], 1,2,1 );
imageplot( clamp(Mwav), ['Wavelets, SNR=' num2str(ewav,3) 'dB'], 1,2,2 );
        
disp(['Delta=' num2str(egeod-ewav) 'dB']);


options.base_str = [rep name '-'];
save_image(clamp({M M1T Mwav}), {'original' ['geod-m' num2str(m) '-snr' num2str(round(10*egeod))] ...
    ['wav-m' num2str(m) '-snr' num2str(round(10*ewav))]}, options);


clf;
plot_triangulation(vertex,faces, M, options);
saveas(gcf, [rep name '-triangulation.eps'], 'eps');