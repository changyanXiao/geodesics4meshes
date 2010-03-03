function plot2D(Im, options)
%--------------------------------------------------------------------------
%   show a 2D image with options for results visu
%   
%   options : - options.fignum      : figure number
%             - options.zoom        : display zoom percent
%             - options.colormap    : colormap
%             - options.colorbar    : horizontal or vertical
%             - options.levelSets   : if so, number of level sets
%             - options.min         : min value for Im to view
%             - options.max         : max value for Im to view
%--------------------------------------------------------------------------
options.null = 0;

if isfield(options, 'fignum')
    fignum = options.fignum;
else
    fignum = 1;
end;

if isfield(options, 'zoom')
    zoom = options.zoom;
else
    zoom = 100;
end;

if isfield(options, 'colormap')
    colmap = options.colormap;
else
    colmap = jet;
end;

if isfield(options, 'colorbar')
    colbar = 'vertical';
else
    colbar = 'horizontal';
end;

if isfield(options, 'min')
    MIN = options.min;
else
    MIN = mmin(Im);
end;

if isfield(options, 'max')
    MAX = options.max;
else
    MAX = mmax(Im);
end;
%--------------------------------------------------------------------------
figure(fignum);
set(gcf,'color',[1 1 1]);
imshow(Im,  [MIN MAX], 'InitialMagnification', zoom);
colormap(colmap);
if isfield(options, 'levelSets')
    hold on;
    contour(Im, options.levelSets, 'Color', 'm');
end;
colorbar(colbar);
axis on;
set(gca,'XTickLabel',{''})
set(gca,'XTick',[]);
set(gca,'YTickLabel',{''})
set(gca,'YTick',[]);
set(gca,'FontName','Times');
    
