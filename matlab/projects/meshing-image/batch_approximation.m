displist = [];


name = 'boat';
n = 256;
% nbr vertex    
m = 600;

metric = 'hessian';
metric = 'structure';

use_lloyd = 1;
epsilon = 1;
ntests = 4;
alpha_list = linspace(1.5,1.5,ntests);
sigma_list = linspace(2,12,ntests)*n/512;
err = [];

rep = ['../../results/image-approx/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

faces_svg = {};
vertex_svg = {};
for ialpha=1:ntests
    alpha = alpha_list(ialpha);
    sigma_structure = sigma_list(ialpha);
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

str_add = '';
clf;
if std(alpha_list)>0
    h = plot(alpha_list,err, '.-');
    title('alpha variation');
    str_add = 'alpha';
else
    h = plot(alpha_list,err, '.-');
    title('sigma variation');
    str_add = 'sigma';
end
set(h, 'LineWidth', 2);
axis tight;
saveas(gcf, [rep name '-m' num2str(m) '-' metric '-lloyd' num2str(use_lloyd) '-' str_add '.png'], 'png');

%% 
% Record optimal triangulation.

[tmp,i] = max(err);
write_off([rep name '-m' num2str(m) '-' metric '-lloyd' num2str(use_lloyd) '.off'], [vertex_svg{i}; zeros(1,m)], faces_svg{i});
