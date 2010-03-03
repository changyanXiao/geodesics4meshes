function [vertex,face] = read_ov(filename)

% read_ov - read data from OV file.
%
%   [vertex,face] = read_ov(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   OV format is "Object Table" format that stores corners. See
%       3D Compression Made Simple: Edgebreaker on a Corner-Table 
%       Jarek Rossignac, Alla Safonova, Andrzej Szymczak  
%   for more invormations.
%
%   Copyright (c) 2008 Gabriel Peyré

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);
[a,str] = strtok(str); nface = str2num(a);

%% retrieve corners %%
[C,cnt] = fscanf(fid,'%f %f', 2*3*nface);
if cnt~=2*3*nface
    warning('Problem in reading corners.');
end
C = reshape(C, [2 3 nface])+1;

%% retieve vertices %%

nvert = [];
while isempty(nvert)
    str = fgets(fid);
    [a,str] = strtok(str); nvert = str2num(a);
end

[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
vertex = reshape(A, 3, cnt/3);

fclose(fid);

%% retrieve face list %%

face = squeeze(C(1,:,:));

