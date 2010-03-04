%PRINT2IM  Prints cropped, anti-aliased figures to image file and/or array
%
% Examples:
%   print2im filename
%   print2im(..., fig_handle)
%   print2im(..., '-alpha')
%   A = print2im(...)
%   [A alpha] = print2im(...)
%
% This function saves and/or outputs a figure as an image nicely, without
% the need to specify multiple options. The image is output at the
% resolution of the figure as it appears on the screen - get your figure
% exactly as you want it, and print!
%
% This function is perfect for saving 3-d graphics for publishing in
% scientific papers. It improves on MATLAB's print command (using default
% options) in these ways:
%   - The figure borders are cropped
%   - Graphics are anti-aliased
%   - Background transparency can be saved (with the png file format) or
%     output as a greyscale alphamatte.
%
% Note that it is preferable to save 2-d graphics as pdfs (using PRINT_PDF)
% if there are vector graphics (e.g. lines, axes) involved, or using the
% original image data (output using SC) if there are only bitmaps (i.e.
% image data without axes) involved, in order to avoid reducing the quality
% of any 2-d bitmaps.
%
%IN:
%   filename - string containing the name (not path) of the file the figure
%              is to be saved as, including the extension, which indicates
%              the file format the image will be saved as. If there is no
%              valid extension (".png", ".tif", ".jpg" or ".bmp") then a
%              ".png" extension is added. The figure is saved in the
%              current directory.
%   fig_handle - The handle of the figure to be saved. Default: current
%                figure.
%   '-alpha' - option specifying that the figure background should be
%              transparent. This option requires that the output file
%              format is png; other formats are changed to png. Note that
%              this option does not automatically set axes background
%              colours to transparent. This can be done using:
%                  set(gca, 'Color', 'none');
%
%OUT:
%   A - MxNxC uint8 output image. C is 1 if the figure is greyscale,
%       otherwise it is 3.
%   alpha - MxN alphamatte for the figure, in the range [0 1]. This assumes
%           the '-alpha' option is given (even if it isn't).
%
% Copyright (C) Oliver Woodford 2008-2009

% This function is inspired by Anders Brun's MYAA (fex id: 20979) and my
% own PRINT_PDF (fex id: 22018).

% $Id: print2im.m,v 1.10 2009/02/12 14:32:25 ojw Exp $

function [A alpha] = print2im(varargin)
% Parse the inputs
[name fig alpha] = parse_inputs(varargin, nargout);
% Is a transparent background required?
if alpha
    % Set the background colour to a colour that hopefully nothing else is
    old_vals = get(fig, {'Color', 'InvertHardcopy'});
    tcol = 255 - [13 3 17];
    set(fig, 'Color', tcol / 255, 'InvertHardcopy', 'off');
    % Print large version to array
    A = print2array(fig, 4);
    % Set the background colour back to normal
    set(fig, 'Color', old_vals{1}, 'InvertHardcopy', old_vals{2});
    % Extract transparent pixels and crop the background
    [A alpha] = crop_background(A, tcol);
    % Set background pixels which will have non-zero alpha to the nearest
    % foreground colour, and paint others white
    A = inpaint_background(A, alpha);
    % Downscale the alphamatte
    alpha = quarter_size(single(alpha), 0);
else
    % Print large version to array
    A = print2array(fig, 4);
    % Crop the background
    A = crop_background(A, [255 255 255]);
    alpha = [];
end
% Downscale
A = quarter_size(A, 255);
% Check if the image is greyscale
if size(A, 3) == 3 && ...
   all(reshape(A(:,:,1) == A(:,:,2), [], 1)) && ...
   all(reshape(A(:,:,2) == A(:,:,3), [], 1))
    A = A(:,:,1); % Save only one channel for 8-bit output
end
if ~isempty(name)
    % Construct the filename
    if numel(name) < 5 || ~any(strcmpi(name(end-3:end), {'.png', '.tif', '.jpg', '.bmp'}))
        name = [name '.png']; % Add the default extension
    end
    name = [cd filesep name]; % Add the path to the current directory
    % Save to file
    if isempty(alpha)
        imwrite(A, name);
    else
        if ~strcmpi(name(end-2:end), 'png')
            name(end-2:end) = 'png';
            warning([mfilename ':IncompatibleExtension'], 'File type (and extension) changed to ''png'' to support transparency');
        end
        imwrite(A, name, 'Alpha', alpha);
    end
end
if nargout < 1
    % Don't print A accidentally
    clear A
end
return

function [name fig alpha] = parse_inputs(inputs, nout)
% Set the defaults
name = [];
fig = gcf;
alpha = false;
% Go through the input arguments
for a = 1:numel(inputs)
    if isnumeric(inputs{a}) && ~isempty(inputs{a})
        fig = inputs{a}(1);
    elseif ischar(inputs{a})
        if inputs{a}(1) == '-'
            switch lower(inputs{a})
                case '-alpha'
                    alpha = true;
                otherwise
                    error('Input argument ''%s'' not recognised', inputs{a});
            end
        else
            name = inputs{a};
        end
    end
end
% Need a filename if no output desired
if isempty(name) && nout < 1
    error('No filename provided')
end
% Need to generate an alphamatte if it is expected as an output
alpha = alpha | nout > 1;
return

function A = print2array(fig, res)
% Generate default input arguments, if needed
if nargin < 2
    res = 1;
    if nargin < 1
        fig = gcf;
    end
end
% Set paper size
set(fig, 'PaperPositionMode', 'auto');
% Generate temporary file name
tmp_nam = [tempname '.tif'];
% Print to tiff file
print(fig, '-opengl', ['-r' num2str(get(0, 'ScreenPixelsPerInch')*res)], '-dtiff', tmp_nam);
% Read in the printed file
A = imread(tmp_nam);
% Delete the file
delete(tmp_nam);
return

function A = quarter_size(A, padval)
% Downsample an image by a factor of 4
try
    % Faster, but requires image processing toolbox
    A = imresize(A, 1/4, 'bilinear');
catch
    % No image processing toolbox - resize manually
    % Lowpass filter - use Gaussian (sigma: 1.7) as is separable, so faster
    filt = single([0.0148395 0.0498173 0.118323 0.198829 0.236384 0.198829 0.118323 0.0498173 0.0148395]);
    B = repmat(single(padval), [size(A, 1) size(A, 2)] + 8);
    for a = 1:size(A, 3)
        B(5:end-4,5:end-4) = A(:,:,a);
        A(:,:,a) = conv2(filt, filt', B, 'valid');
    end
    clear B
    % Subsample
    A = A(2:4:end,2:4:end,:);
end
return

function [A alpha] = crop_background(A, bcol)
% Map the foreground pixels
alpha = A(:,:,1) ~= bcol(1) | A(:,:,2) ~= bcol(2) | A(:,:,3) ~= bcol(3);
% Crop the background
N = any(alpha, 1);
M = any(alpha, 2);
M = find(M, 1):find(M, 1, 'last');
N = find(N, 1):find(N, 1, 'last');
A = A(M,N,:);
if nargout > 1
    % Crop the map
    alpha = alpha(M,N);
end
return

function A = inpaint_background(A, alpha)
% Inpaint some of the background pixels with the colour of the nearest
% foreground neighbour
% Create neighbourhood
[Y X] = ndgrid(-4:4, -4:4);
X = Y .^ 2 + X .^ 2;
[X I] = sort(X(:));
X(I) = 2 .^ (numel(I):-1:1); % Use powers of 2
X = reshape(single(X), 9, 9);
X = X(end:-1:1,end:-1:1); % Flip for convolution
% Convolve with the mask & compute closest neighbour
M = conv2(single(alpha), X, 'same');
J = find(M ~= 0 & ~alpha);
[M M] = log2(M(J));
% Compute the index of the closest neighbour
[Y X] = ndgrid(-4:4, (-4:4)*size(alpha, 1));
X = X + Y;
X = X(I);
M = X(numel(X) + 2 - M) + J;
% Reshape for colour transfer
sz = size(A);
A = reshape(A, [sz(1)*sz(2) sz(3)]);
% Set background pixels to white (in case figure is greyscale)
A(~alpha,:) = 255;
% Change background colour to closest foreground colour
A(J,:) = A(M,:);
% Reshape back
A = reshape(A, sz);
return