function viewBySlice(IM)
%
%
%
figure;
[nx, ny, nz] = size(IM);
m = mmin(IM); M = mmax(IM);
for k =1:nz
    imshow(IM(:,:,k), [m M]);
    title(int2str(k));
    pause;
end;