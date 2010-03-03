%PIX2ANIGIF Combines sequence of numbered pictures into an animated GIF
%
%   pix2anigif(inputFileRoot,inputFileType,firstImage,lastImage,...
%          outputFileRoot,numberPattern,framesPerSecond,loopCount,...
%          displayWhileCreating)
%
%   - inputFileRoot: String containing the path and root file name for
%                    the sequence of input images.  For example, if the
%                    files are of the form 'frames\input_001.png', this
%                    variable should be 'frames\input_'
%   - inputFileType: String containing the extension of the file type.
%                    For example, 'bmp', 'png', 'gif', 'jpg'.  Any image
%                    supported by imread will work.
%   - firstImage: Index of the first image to use.  Typically, 0 or 1.
%   - lastImage: Index of the last image to use.  It must be greater than
%                firstImage.
%   - numberPattern: (optional) Format string for indices.  For instance,
%                    for indices 000, 001, 002, ... numberString should be
%                    '%03d' where zero specifies the zero padding, 3
%                    specifies the width, and d is the integer.  If the
%                    numbers are 0, 1, 2, ... numberstring should be '%d'
%                    so that extra zeros and spaces aren't inserted.
%                    Default is '%03d'.
%   - outputFileRoot: (optional) String containing the output file name,
%                     without the extension.  For example,
%                     'animatedOutput'.  Directories can be included in
%                     the string.  Default is 'out'.
%                *** NOTE: The output file is automatically overwritten.
%   - framesPerSecond: (optional) Desired frame rate for playback.  Note,
%                      the top speed is limited by the playback software.
%                      Default is 10 fps.
%   - loopCount: (optional) Number of times to loop during playback.  To
%                loop continuously (default) set this value to Inf.
%   - displayWhileCreating: (optional) This is a boolean input.  When it 
%                           is true, images are displayed as they are
%                           processed.  Default is false.
%
%   Automatically creates animated GIFs from a numbered sequence of
%   images using Matlab's imwrite function.

%   (c) Gabe Hoffmann, gabe.hoffmann@gmail.com
%   Written 8/31/2007

function pix2anigif(inputFileRoot,inputFileType,firstImage,lastImage,...
                outputFileRoot,numberPattern,framesPerSecond,loopCount,...
                displayWhileCreating)

% Verify correct number of arguments
error(nargchk(4,9,nargin));

if firstImage >= lastImage
    error('lastImage must be greater than firstImage');
end

% Set default values for optional arguments
if nargin < 5
    outputFileRoot = 'out';
end
if nargin < 6
    numberPattern = '%03d';
end
if nargin < 7
    framesPerSecond = 10;
end
if nargin < 8
    loopCount = Inf;
end
if nargin < 9
    displayWhileCreating = false;
end
for i = firstImage:lastImage
    % Load next frame
    [icdata,imap] = imread([inputFileRoot,sprintf(numberPattern,i),...
                            '.',inputFileType]);
    imap = jet(256);
	% Add frames to GIF.  For the first image, overwrite any existing
    % images
    if i==firstImage
        imwrite(icdata,imap,[outputFileRoot,'.gif'],'gif',...
            'DelayTime',1/framesPerSecond, 'LoopCount',loopCount,...
            'WriteMode','overwrite')
    else
        imwrite(icdata,imap,[outputFileRoot,'.gif'],'gif','WriteMode','append',...
            'DelayTime',1/framesPerSecond)
    end
    
    % Show images in figure if displayWhileCreating == true
    if displayWhileCreating
        imshow(icdata,imap)
        drawnow
    end
end

% Display completion
fprintf('Finished creating %s\n',[outputFileRoot,'.gif']);