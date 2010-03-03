% test for 3x3 determinant

p = 100000;
A = randn(3,3,p);

tic;
d = det3(A);
toc;

tic;
d1 = zeros(p,1);
for i=1:p
    d1(i) = det(A(:,:,i));
end
toc;

norm(d-d1)