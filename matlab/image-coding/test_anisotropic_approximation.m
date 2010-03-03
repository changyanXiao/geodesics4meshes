%% Test for m-triangle approximation

test_generic = 1;
test_alpha = 0;
test_m = 0;

rep = 'results/anisotropic-approximation/';
if not(exist(rep))
    mkdir(rep);
end

%% Load the image
name = 'cameraman';
n = 256;
M = load_image(name,n);
M = rescale(crop(M,n));

%% Parameter for the tensor field remapping
epsilon = 1e-6;
%% Parameter for interpolation
options.remove_nan = 1;
%% Structure tensor extraction
sigma1 = 1.5;   % for gradient computation
sigma2 = 2;     % for filtering the tensor field
options.use_anisotropic = 0;
T0 = compute_structure_tensor(M, sigma1, sigma2);
[e1,e2,l1,l2] = perform_tensor_decomp(T0);

% U = perform_tensor_mapping(T,+1);

%% Generic display
if test_generic
    epslist = [1e-8 1e-6 1e-4 1e-3];
    for iepsilon = 1:length(epslist)
    for alpha = [1]
        for m=[600]
            epsilon = epslist(iepsilon);
            
            str = ['-m' num2str(m) '-alpha' num2str(10*alpha) '-eps' num2str(iepsilon) ];
            
            % Structure tensor remapping
            % l1 > l2, e1 is alligned with contour
            T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);
            
            % Farthest seeding
            [U, dUx,dUy,V, vertex,faces] = FPSCAniso2D([1;1], T, m, 0);
            faces = faces+1;
            
            options.edgewidth = 1.5;
            options.vertexsize = 15;
            
            % display
            clf;
            plot_triangulation(vertex,faces, M, options);
            % saveas(gcf, [rep name str '.eps'], 'eps');
            saveas(gcf, [rep name str '-triang.png'], 'png');
            
            % approximated image
            options.verb = 0;
            v1 = compute_orthoproj_triangulation(vertex, faces, M);
            M1 = griddata_arbitrary(faces,vertex,v1,n, options);
            
            % save image
            options.base_str = [rep name];
            % save_image(clamp({M M1}), {'-original' str}, options);
            % SNR display
            disp(['m=' num2str(m) ' alpha=' num2str(alpha) ...
                ' epsilon=' num2str(epsilon) ' SNR=' num2str(snr(M,M1)) ]);
        end
    end
    end    
end
%%%%%%


%% testing the influence of anisotropy
if test_alpha
    ntests = 14;
    m = 400;
    alpha_list = linspace(.5,5,ntests);
    err = [];
    for i=1:ntests
        alpha = alpha_list(i);
        T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);
        [U, dUx,dUy,V, vertex,faces] = FPSCAniso2D([1;1], T, m, 0);
        faces = faces+1;
        % approximate
        v1 = compute_orthoproj_triangulation(vertex, faces, M);
        options.remove_nan = 1;
        M1 = griddata_arbitrary(faces,vertex,v1,n, options);
        % keep error
        err(i) = snr(M, M1);
        % save_image(clamp({M M1}), {'original' ['m-' num2str(m)]}, options);
    end
    
    clf;
    h = plot(alpha_list, err, 'k');
    axis tight; box on;
    set(h, 'LineWidth', 2);
    saveas(gcf, [rep name '-influ-alpha.eps'], 'eps');
end
%%%%%%%


%% testing the error decay %%
if test_m
    
    ntests = 6;
    mlist = linspace(100, 1000, ntests);
    
    alphalist = linspace(.4, 1, 4);
    err_geod = []; mlist_real = [];
    err_eucl = [];
    for i=1:ntests
        m = mlist(i);
        disp(['--- m=' num2str(m) ' ----']);
        % approximate with anisotropic geodesic delaunay
        for k=1:length(alphalist)
            alpha = alphalist(k);
            T = perform_tensor_recomp(e2, e1, (l1+epsilon).^alpha, (l2+epsilon).^alpha);
            [U, dUx,dUy,V, vertex,faces] = FPSCAniso2D([1;1], T, m, 0); faces = faces+1;
            mlist_real(i,k) = size(vertex,2);
            % geodesic
            v1 = compute_orthoproj_triangulation(vertex, faces, M);
            M1 = griddata_arbitrary(faces,vertex,v1,n, options);
            err_geod(i,k) = snr(M, M1);
            % euclidean
            faces_euc = compute_delaunay(vertex);
            v1 = compute_orthoproj_triangulation(vertex, faces_euc, M);
            M1 = griddata_arbitrary(faces_euc,vertex,v1,n, options);
            err_eucl(i,k) = snr(M, M1);
        end
    end
    klist = [1 3 4];
    gr = {'' '.' '*' 'o' '^' 'v'};
    ms = [1 20 10 10 10 10]; lw = 2;
    lgd= {};
    clf;
    hold on;
    for i=1:length(klist)
        k = klist(i);
        lgd{i} = ['\alpha=' num2str(alphalist(k))];
        h = plot( mlist_real(:,k), err_geod(:,k), ['k' gr{i} '-']);
        set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms(i));
    end
    for i=1:length(klist)
        k = klist(i);
        plot( mlist_real(:,k), err_eucl(:,k), ['k' gr{i} '--']);
        set(h, 'LineWidth', 1); set(h, 'MarkerSize', ms(i));
    end
    legend(lgd);
    axis tight; box on;
    saveas(gcf, [rep name '-influ-m.eps'], 'eps');    
end
%%%%%
