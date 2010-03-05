% test for coding of a sparse set of points

n = 100;
k=(0:n)';

e1 = k.*log2(n);
e2 = k.*log2(n./k)+(n-k).*log2(n./(n-k));
e2(end) = 0;
e1(1)= 0; e2(1)= 0;

e3 = [];
for i=1:length(k)
    e3(i) = log2(nchoosek(n,k(i)));
end
e3 = e3(:);

clf;
plot(k/n, [e1 e2 e3]);
axis tight;
legend('independant', 'shanon', 'nchoosek');