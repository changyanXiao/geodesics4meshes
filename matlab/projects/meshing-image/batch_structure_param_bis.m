add_base_paths;

default name = ['piecewise-quadratic'];

default n = 256;
% nbr vertex    
default m = 800;

metric = 'hessian';
metric = 'structure';
displist = [];

ntests = 8; 

sigma_structure = 12*n/512;
alpha_list = linspace(2,4,ntests);
epsilon_list = linspace(1e-5,1e-2,ntests);


use_lloyd = 1;
err = zeros(ntests,ntests);

rep = ['../../results/image-approx/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

faces_svg = {};
vertex_svg = {};
for ialpha=1:ntests
for iepsilon=1:ntests
    alpha = alpha_list(ialpha);
    epsilon = epsilon_list(iepsilon);
    test_approximation;
    err(ialpha,iepsilon) = snr(M,M1);
    %%
    % Display
    if 0
        options.col = 'b-'; options.ms = 0; options.ps = 1;
        clf; hold on;
        imageplot(M);        
        plot_graph(triangulation2adjacency(faces),vertex(2:-1:1,:), options);
        h = plot([1 n n 1 1], [1 1 n n 1], 'b.-');
        set(h, 'LineWidth', 1);
        axis ij;
        saveas(gcf, [rep name '-m' num2str num2str(alpha)]);
    end
    faces_svg{ialpha,iepsilon} = faces;
    vertex_svg{ialpha,iepsilon} = vertex;
end
end

return;
%%
% Display

switch test_type
    case 'alpha'
        x = alpha_list;
        tlt = ['epsilon=' num2str(epsilon_list(1)) ', sigma=' num2str(sigma_list(1)) ];
    case 'sigma'
        x = sigma_list;
        tlt = ['epsilon=' num2str(epsilon_list(1)) ', alpha=' num2str(alpha_list(1)) ];
    case 'epsilon'
        x = epsilon_list;
        tlt = ['alpha=' num2str(alpha_list(1)) ', sigma=' num2str(sigma_list(1)) ];
end

clf;
h = plot(x,err, '.-');
title(tlt);
xlabel(test_type);
set(h, 'LineWidth', 2);
axis tight;
saveas(gcf, [rep name '-m' num2str(m) '-' metric '-lloyd' num2str(use_lloyd) '-' test_type '.png'], 'png');

return;

%% 
% Record optimal triangulation.

[tmp,i] = max(err);
write_off([rep name '-m' num2str(m) '-' metric '-lloyd' num2str(use_lloyd) '.off'], [vertex_svg{i}; zeros(1,m)], faces_svg{i});
