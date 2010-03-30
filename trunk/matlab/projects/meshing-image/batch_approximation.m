%% batch for several images, using well chosen parameters (alpha,sigma)

name_list = {'boat' 'game' 'piecewise-quadratic'};
% name_list = {'piecewise-quadratic'};


rep = ['../../results/image-approx/svg/'];
if not(exist(rep))
    mkdir(rep);
end

% - Boat avec 7000 sommets (PSNR: 31.83 dB)
% m=2000 SNR=20.6 pour alpha=1.6, epsilon=1e-3, sigma=6*n/512
% - Game avec 6000 sommets (PSNR: 36.54 dB)

% - PiecewiseQuadratic avec 800 sommets (PSNR: 42.85 dB)
% SNR 30 env, pour alpha=3, epsilon=1e-3, sigma=4*n/512

n = 256*2;
metric = 'structure';
use_lloyd = 0;
lloyd_refresh = 100;
displist = [];


mlist = round([7000 6000 800]);
alpha_list = [1.6 1.6 3];
sigma_list = [6 6 4]*n/512;
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