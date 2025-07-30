%%TEST

v = [1 2 3];
n = length(v);
a = zeros(1, n);
for k = n:-1:1
    a(k) = v(n-k+1);
end

a = a(1:end-1);

Phi = [1 2; 3 4; 5 6; 7 8];

R = Phi(end-2:end, :);