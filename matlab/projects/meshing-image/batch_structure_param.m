add_base_paths;

default name = ['boat'];
% name = 'piecewise-quadratic';

default n = 256;
% nbr vertex    
default m = 800;

metric = 'hessian';
metric = 'structure';
displist = [];

ntests = 6; 

alpha_default = 1.6;
epsilon_default = .1*1e-2;
sigma_default = 5*n/512;

% default values
switch name
    case 'boat'
        alpha_default = 1.6;
        epsilon_default = 1e-3;
        sigma_default = 8*n/512;
        alpha_list = linspace(1.4,2.2,ntests);
        epsilon_list = linspace(1e-5,1e-3,ntests);
        sigma_list = linspace(4,20,ntests)*n/512;           
    case 'piecewise-quadratic'
        % SNR 27 env
        alpha_default = 2;
        epsilon_default = .001;
        sigma_default = 12*n/512; 
        alpha_list = linspace(1.5,3.5,ntests);
        epsilon_list = linspace(1e-5,1e-3,ntests);
        sigma_list = linspace(4,20,ntests)*n/512;   
end


default test_type = ['alpha'];

switch test_type
    case 'alpha'
        epsilon_list = ones(ntests,1)*epsilon_default;
        sigma_list = ones(ntests,1)*sigma_default;
    case 'epsilon'
        alpha_list = ones(ntests,1)*alpha_default;
        sigma_list = ones(ntests,1)*sigma_default;
    case 'sigma'
        epsilon_list = ones(ntests,1)*epsilon_default;
        alpha_list = ones(ntests,1)*alpha_default;
end

use_lloyd = 0;
err = [];

rep = ['../../results/image-approx/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

faces_svg = {};
vertex_svg = {};
for itest=1:ntests
    alpha = alpha_list(itest);
    sigma_structure = sigma_list(itest);
    epsilon = epsilon_list(itest);
    test_approximation;
    err(end+1) = snr(M,M1);
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
    faces_svg{end+1} = faces;
    vertex_svg{end+1} = vertex;
end

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
