function save_gcf(filename, save_eps,save_png)

% save_gcf - save gcf in png and/or eps
%
%   save_gcf(filename, save_eps,save_png);
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<2
    save_eps = 1;
end
if nargin<3
    save_png = 1;
end

if save_eps
    saveas(gcf, [filename '.eps'], 'epsc');
end
if save_png
    saveas(gcf, [filename '.png'], 'png');
end