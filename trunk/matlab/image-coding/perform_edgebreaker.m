function [face, vertex] = perform_edgebreaker(face, vertex, direction, options)

% perform_edgebreaker - wrapper to edgebreaker topology compression
%
% Compression of topology + geometry:
%   [stream,vertexC] = perform_edgebreaker(face, vertex, +1, options);
% Compression of topology alone:
%   stream = perform_edgebreaker(face, [], +1, options);
%
% Decompression:
%   [face, vertex] = perform_edgebreaker(stream, vertexC, -1, options);
% Compression of topology alone:
%   face = perform_edgebreaker(stream, [], -1, options);
%
%   IMPORTANT: if your surface is not a close manifold, then you must set
%   options.patch = 1. In this case, only one boundary is allowed, the
%   others (smaller) holes are filled.
%
%   Set options.coding='binary' to convert the integer stream to a binary stream.
%
%   This is wrapper to the Edgebreaker code that can be downloaded from
%   here
%       http://www.gvu.gatech.edu/~jarek/edgebreaker/eb/EBsoftware.html
%
%   The executables EBCompression and EBDecompression needs to be
%   available.
%   The location of the EB executables can be given in options.EBdir
%   (default './').
%
%   Copyright (c) 2008 Gabriel Peyre


options.null = 0;
binary = getoptions(options, 'binary', 0);
EBdir = getoptions(options, 'EBdir', './');
patch = getoptions(options, 'patch', 0);
coding = getoptions(options, 'coding', 'int');

format = 'ASKII';
if binary==1
    format = 'BINARY';
end

surf_type = 'MANIFOLD';
if patch==1
    surf_type = 'TPATCH';
end

ov_filename = 'tmp';
ov_filename = [ov_filename '.ov'];
output_dir = '.';
clers_file = 'clers.txt';
handles_file = 'handles.txt';
vertices_file = 'vertices.txt';

StartCorner = 0;

% ov_filename = 'OVTable_2Handle.ov';

%   StartCorner: Is the corner where to start EBCompression.
% -           If  MeshType is TPATCH it should be a corner corresponding to the dummy vertex.
% -           If  MeshType is MANIFOLD, it can be any corner, but since the triangles incident
% on  StartCorner are not stored it is advantageous to pass a corner with the maximum number of incident triangles.

if direction==1
    if isempty(vertex)
        nvert = max(face(:));
        vertex = zeros(3, nvert);
    end
    nvert = size(vertex,2);
    if size(vertex,1)==2
        vertex = cat(1, vertex, zeros(1,nvert));
    end
    dummy = write_ov(ov_filename, vertex, face);
    if dummy>0
        StartCorner = dummy-1;
    end
    %% CODING %%
    system([EBdir 'EBCompression ' ov_filename ' ' surf_type ' ' num2str(StartCorner)  ' ' format ' > tmp']);
    system(['rm ' ov_filename]);
    system(['rm tmp']);
    % recover the information from clers.txt
    fid = fopen(clers_file, 'rt');  
    if fid<=0
        error('Problem with EB');
    end  
    [t,cnt] = fscanf(fid,'%s', 1);
    [a,cnt] = fscanf(fid,'%d', 1);
    [C,cnt] = fscanf(fid,'%s', Inf);
    % recover the information from vertices.txt
    fid = fopen(vertices_file, 'rt');  
    if fid<=0
        error('Problem with EB');
    end
    [vertex,cnt] = fscanf(fid,'%f %f %f', Inf);
    vertex = reshape(vertex, [3 length(vertex)/3]);
    % store as integers
    stream = clers2num(C);
    % store as bits
    if strcmp(coding, 'binary')
        stream = coding_int2bin(stream);
    end
    % add integers
    stream = [a; stream];    
    %
    fclose(fid);
    system(['rm ' clers_file]);
    system(['rm ' handles_file]);
    system(['rm ' vertices_file]);
    % output
    face = stream;
else
    stream = face;
    a = stream(1);
    % b = stream(2);
    stream = stream(2:end);
    if strcmp(coding, 'binary')
        stream = coding_int2bin(stream);
    end
    C = num2clers(stream);  
    % write the information to clers.txt
    fid = fopen(clers_file, 'wt');   
    if fid<=0
        error('Problem with EB');
    end   
    fprintf(fid, '%s\n', surf_type);
    fprintf(fid, '%d\n', a);  
    for k=1:length(C)
        fprintf(fid, '%c\n', C(k));
    end
    fclose(fid);
    % write the information to vertex
    fid = fopen(vertices_file, 'wt');  
    if fid<=0
        error('Problem with EB');
    end
    fprintf(fid, '%f %f %f\n', vertex);   
    fclose(fid);
    % write the information to handles : nothing
    fid = fopen(handles_file, 'wt');  
    fclose(fid);    
    % run EB
    system([EBdir 'EBDecompression ' output_dir ' ' ov_filename ' ' format ' > tmp']);
    % recover data from the file
    [vertex,face] = read_ov(ov_filename);
    % remove tmp files
    system(['rm ' ov_filename]);
    system(['rm tmp']);
    % remove dummy vertex, whose Id is 1
    if patch==1
        dummy = size(vertex,1);        
        I = find(sum(face==1)==0); face = face(:,I)-1;
        vertex = vertex(:,2:end);
    end
    if sum(abs(vertex(3,:)))<1e-9
        vertex = vertex(1:2,:);
    end
end

%%
function stream = clers2num(C)

stream = zeros(length(C), 1);
for i=1:length(C)
    switch C(i)
        case 'C'
            c = 0;
        case 'L'
            c = 1;
        case 'E'
            c = 2;
        case 'R'
            c = 3;
        case 'S'
            c = 4;
        otherwise
            error('Unknown symbol');
    end
    stream(i) = c;
end

%%
function C = num2clers(stream)

C = [];
for i=1:length(stream)
    switch stream(i)
        case 0
            c = 'C';
        case 1
            c = 'L';
        case 2
            c = 'E';
        case 3
            c = 'R';
        case 4
            c = 'S';
        otherwise
            error('Unknown symbol');
    end
    C = [C c];
end