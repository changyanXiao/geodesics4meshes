function dummy = write_ov(filename, vertex, face);

% write_off - write a mesh to a OV file
%
%   dummy = write_ov(filename, vertex, face);
%
%   vertex must be of size [n,3]
%   face must be of size [p,3]
%
%   OV format is "Object Table" format that stores corners. See
%       3D Compression Made Simple: Edgebreaker on a Corner-Table 
%       Jarek Rossignac, Alla Safonova, Andrzej Szymczak  
%   for more invormations.
%
%   Copyright (c) 2008 Gabriel Peyré


if size(vertex,1)~=3
    vertex=vertex';
end
if size(vertex,1)~=3
    error('vertex does not have the correct format.');
end

if size(face,1)~=3
    face=face';
end
if size(face,1)~=3
    error('face does not have the correct format.');
end

nface = size(face,2);
nvert = size(vertex,2);

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

% compute boundary in the mesh and add dummy vertices
options.verb = 0;
bound = compute_boundary(face, options);
dummy = -1;
if not(isempty(bound))
    boundary = bound(:)';
    nbound = length(bound);
    face(:,end+1:end+nbound) = [bound; bound([2:end 1]); ones(1,nbound)*(nvert+1)];
    nface = nface + nbound;
    % add dummy vertex
    nvert = nvert+1;
    vertex(:,end+1) = [0;0;0];
    % corner corresponding to dummy vertex
    dummy = nface*3;
end

% find the opposite point for each face elements
options.verb = 0;
e2f = compute_edge_face_ring(face);
faceo = zeros(3,nface);
for f=1:nface
    for k=1:3
        i = face(k,f);
        i1 = face( mod(k,3)+1  ,f);
        i2 = face( mod(k+1,3)+1,f);
        % find oposite face
        f1 = e2f(i1,i2); f2 = e2f(i2,i1);
        if f1<=0 && f2<=0
            error('Problem with mesh');
        end
        if f1==f
            f1 = f2;
        end
        % find number
        k1 = find( face(:,f1)==i1 );
        k2 = find( face(:,f1)==i2 );
        kk = [1 2 3]; kk([k1 k2]) = [];
        % asign
        faceo(k,f) = 3*(f1-1)+kk;        
    end
end

% write number of faces
fprintf(fid, '%d\n', nface );

% write corner list
C = [face(:) faceo(:)]';
fprintf(fid, '%d %d\n', C(:)-1);

% write number of vertices
fprintf(fid, '%d\n', nvert );

% write the vertices
fprintf(fid, '%f %f %f\n', vertex);

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i1,i2] = find_other(f,i)
f(f==i)=[];
if length(f)~=2
    error('Problem');
end
i1 = f(1); i2 = f(2);