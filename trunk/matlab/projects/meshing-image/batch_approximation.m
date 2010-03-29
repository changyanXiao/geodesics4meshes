%% batch for several images, using well chosen parameters (alpha,sigma)

name_list = {'boat' 'game' 'piecewise-quadratic'};
% name_list = {'piecewise-quadratic'};


rep = ['../../results/image-approx/svg/'];
if not(exist(rep))
    mkdir(rep);
end

% - Boat avec 7000 sommets (PSNR: 31.83 dB)
% - Game avec 6000 sommets (PSNR: 36.54 dB)
% - PiecewiseQuadratic avec 800 sommets (PSNR: 42.85 dB)

n = 256*2;
metric = 'structure';
use_lloyd = 1;
lloyd_refresh = 100;
displist = [];


mlist = round([7000 6000 800]);
alpha_list = [1.6 1.6 3];
sigma_list = [4 4 8]*n/512;
epsilon_list = [1e-3 1e-3 1e-3];

for iname = 1:length(name_list)
    
    
    epsilon = epsilon_list(iname);
    alpha = alpha_list(iname);
    sigma_structure = sigma_list(iname);

    name = name_list{iname};
    m = mlist(iname);
    test_approximation;
    % save image
    warning off;
    imwrite(clamp(M1), [rep name '-approx-m' num2str(m) '.png'], 'png');
    warning on;
    % save triangulation
    write_off([rep name '-m' num2str(m) '-triangles-m' num2str(m)  '.off'], [vertex; zeros(1,m)], faces);
    disp([name ', PSNR:' num2str(psnr(M,M1), 3) 'dB']);
end