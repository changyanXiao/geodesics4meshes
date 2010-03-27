%% batch for several images, using well chosen parameters (alpha,sigma)

name_list = {'boat' 'game' 'piecewise-quadratic'};
mlist = round([7000 6000 800]);


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
lloyd_refresh = 50;
epsilon = 1;
alpha = 1.5;
sigma_structure = 10*n/512;
displist = [];

for iname = 1:length(name_list)
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