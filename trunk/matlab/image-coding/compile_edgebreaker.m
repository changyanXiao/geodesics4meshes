% compile edgebreaker

cd edgebreaker;
!g++ EBCompression.cpp -o EBCompression
!g++ EBDecompression.cpp -o EBDecompression
!mv EBDecompression ../
!mv EBCompression ../
cd ..
% 
% 1.        OVTable: The Corner Table containing the connectivity and
%                   geometry of the mesh in OVTable format 
% 2.        MeshType:  MANIFOLD or TPATCH 
% -           MANIFOLD - a consistently oriented, manifold Triangle mesh, without boundary.                                
% -           TPATCH - a manifold mesh with a single hole that has been 
%                   closed with dummy vertices and triangles 
% 3.        StartCorner: Is the corner where to start EBCompression. 
% -           If  MeshType is TPATCH it should be a corner 
%                   corresponding to the dummy vertex. 
% -           If  MeshType is MANIFOLD, it can be any corner, 
%                   but since the triangles incident on  StartCorner are not 
%                   stored it is advantageous to pass a corner with the maximum 
%                   number of incident triangles. 
% 4.        FileFormat: BINARY or ASCII  