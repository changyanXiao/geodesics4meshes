%% Geometrical Coder for Sparse Points
% We show how a contextual coder can improve for the coding of 
% points with clusters.

%% Generate a 2D Set of Points

%% 
% Size of the image.

n = 256;

%% 
% Number of points.

k = round(n^2*.3);

%%
% Generate a binary image, with or without clusters.

name = 'rand';
name = 'grad';
switch name
    case 'rand'
        % binary image
        A = zeros(n);
        sel = randperm(n^2); sel = sel(1:k);
        A(sel)=1;
    case 'grad'
        M = load_image('lena',n);
        G = grad(M);
        G = sum(G.^2, 3);
        [tmp,sel] = sort(G(:), 'descend');
        A = zeros(n);
        A(sel(1:k)) = 1;
end

%% 
% Binary image display.

clf;
imageplot(A);

%% Naive Coding
% A naive coding consist in treating the black/white point as
% generate from a Bernouilli binary distribution with probability
% |k/n^2| of having 1. Then we use the Shanon bound to compute
% the required number of bits.

%%
% Compute the number of bits 
% for a naive independant coding of the points.

N = n*n;
Enaive = k.*log2(N./k)+(N-k).*log2(N./(N-k));

%% Generate a Conditional Coding Context
% One can improve the coding by using a contextual coding.
% We split the pixels into several ground, according to the number of
% white points in a causal neighborhood. If there are 
% clusters of points, this will be likely to improve the coding rate.

%%
% Very important parameter: size of the context.

s = 4;

%%
% Compute the context neighbor mask

H = zeros(2*s+1);
H(1:s,:) = 1;
H(s+1,1:s) = 1;

%% 
% Compute the context.

C = perform_convolution(A,H);

%%
% Display context.

clf;
imageplot(A, 'Image', 1,2,1);
imageplot(C, 'Context', 1,2,2);

%%
% Compute histograms.

Cmax = max(C(:));
t = 0:Cmax;
h = hist(C(:), t);

%%
% Display histograms of the distribution of context values.

clf; 
bar(h);
axis tight;

%% Perform the Contextual Coding
% We use a small number of class and partition the data according
% to the class.

%%
% Important parameter: number of class for the coder
% (increasing |q| will improve the results, but 
% using a too large value is cheating because here
% we are using an asymptotic Shannon entropy bound).

q = 5;

%%
% Compute the classes for the entropy computation.

Ctxt = {}; Ctxt{1} = 0;
for i=1:q
    Cprev = Ctxt{end}(end);
    Cnew = Cprev + round(Cmax/q);
    if i==q
        Cnew = Cmax;
    end
    Ctxt{end+1} = Cprev+1:Cnew;
end

%%
% Compute the number of bits with this conditional coding.

Aclust = zeros(n);
Eclust = 0;
p = [];
for i=1:length(Ctxt)
    c = Ctxt{i};
    if not(isempty(c))
        I = find( C>=c(1) & C<=c(end) );
        Aclust(I) = i;
        p(i) = min(max(sum(A(I)==1) / length(I), 1e-10), 1-1e-10);
        % entropy
        e = - p(i)*log2(p(i)) - (1-p(i))*log2(1-p(i));
        Eclust = Eclust + length(I) * e;
    end
end

%%
% Display the clustering of the points into classes.

clf;
imageplot(Aclust);
colormap jet;

%% 
% Display how the probability |p(i)| of having 1 increase with the class
% number.

clf;
plot(p, '.-');
axis('tight');

%%
% Display the gain.

disp(['Performance gain: '  num2str(100*(Enaive/Eclust-1)) '%']);

