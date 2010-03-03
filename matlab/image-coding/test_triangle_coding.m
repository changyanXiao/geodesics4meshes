% test for coding of position and values 

% algorithm
% - compute the delaunay
% - compute the position quantization that generates approx. the same error

n = 512;
M = load_image('lena', n);
M = rescale(M);

% quantization of values
qval = 255;
% quantization of position
qpos = n;

% sample random triangles
plist = round(linspace(300,5000,10));
ntests = length(plist);

err_original = [];
err_transformed = [];
R_original = [];
R_transformed = [];
vertex = [0 0 1 1; 0 1 0 1]*(n-1)+1;

for i=1:ntests
    disp(['--> Test #' num2str(i) '/' num2str(ntests) '.']);
    p = plist(i);
    dp = p - size(vertex,2);
    
    % generate quantized positions
    V = floor( rand(2,dp)*(qpos-1)+1 );
    V = V/qpos*n;
    vertex(:,end+1:p) = V;
    face = compute_delaunay(vertex+randn(size(vertex))*1e-9);
    

    %%% Non-transformed approximation %%%
    % compute values
    v = compute_orthoproj_triangulation(vertex, face, M);
    v = clamp(v);
    [tmp, vQ] = quantize(v,1/qval); vQ = clamp(vQ,0,1);
    M1 = griddata_arbitrary(face,vertex,vQ,n);
    err_original(i) = snr(M,M1);
    % number of bits = #bits.position + #bits.values
    R_original(i) = 2*p*log2(qpos) + p*log2(qval);
    
    %%% Transformed approximation %%%
    
    % the graph associated to the triangulation
    A = triangulation2adjacency(face);
    
    % target error on the values
    eta = 1/10; % percent of quantization
    eValTgt =  eta * 1/qval * p;
    eVertTgt = eta * n/qpos * 2*p;
    
    % compute an orthogonal basis for re-transforming the position/values
    if 0
        options.symmetrize=1; options.normalize=0;
        L = compute_mesh_laplacian(vertex,face,'combinatorial',options);
        fprintf('Computing eigenvectors ... ');
        [U,S] = eig(full(L));
        fprintf('done.\n');
        uVert = (U'*vertex')';
        uVal  = (U'*v);
    else
        uVal = perform_haar_graph(v, A, +1, options);
        uVert = zeros(2,p);
        for s=1:2
            uVert(s,:) = perform_haar_graph(vertex(s,:)', A, +1, options)';
        end
    end
    
    
    % search for a quantization that provide a target error
    TlistVert = linspace(1e-5,mmax(abs(uVert))/10,2000);
    TlistVal = linspace(1e-5,mmax(abs(uVal))/10,2000);
    eVert = []; eVal = [];
    for k=1:length(Tlist)        
        [uVertI, uVertQ] = quantize(uVert,TlistVert(k));
        [uValI, uValQ] = quantize(uVal,TlistVal(k));
        eVert(k) = norm( uVertQ-uVert, 'fro' );
        eVal(k) = norm( uValQ-uVal, 'fro' );
    end
    [tmp,k] = min( abs(eVertTgt-eVert) ); TVert = TlistVert(k);
    if k==1 || k==length(Tlist), warning('Problem'); end
    [tmp,k] = min( abs(eValTgt-eVal) ); TVal = TlistVal(k);
    if k==1 || k==length(Tlist), warning('Problem'); end
    % quantize position
	[uVertI, uVertQ] = quantize(uVert,TVert);
	[uValI, uValQ] = quantize(uVal,TVal);
    % un-do the transform
    if 0
        vertexQ = (U*uVertQ')';
        vQ = U*uValQ;
    else
        vQ = perform_haar_graph(uValQ, A, -1, options);
        vertexQ = zeros(2,p);
        for s=1:2
            vertexQ(s,:) = perform_haar_graph(uVertQ(s,:)', A, -1, options)';
        end
    end
    % to check correctness
    % norm(v-vQ, 'fro')/eValTgt
    % norm(vertex-vertexQ, 'fro')/eVertTgt
    % interpolation with approximation position and values
    M2 = griddata_arbitrary(face,vertexQ,v,n);
    err_transformed(i) = snr(M,M2);
    R_transformed(i) = 2*p*compute_entropy(uVertI(:)) + p*compute_entropy(uValI(:));
    
    % real coding - not cheating
    options.coder_type = 1;
    [y,Rvert] = perform_arithmetic_coding(uVertI(:), +1, options);
    [y,Rval] = perform_arithmetic_coding(uValI(:), +1, options);
    R_transformedC(i) = Rvert + Rval;
    
end

clf;
plot(R_original/n^2,err_original, R_transformed/n^2,err_transformed, R_transformedC/n^2,err_transformed);
xlabel('#bits/pixel'); ylabel('SNR');
legend('Original', 'Transformed (entropy)', 'Transformed (arithmetic)');