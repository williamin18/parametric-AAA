function b = multi_r1_times_TT(A,x)
[n_samples,~] = size(A{1});
[n_d,N,R] = TTsizes(x);
W = A{1}*x{1};
for i = 2:n_d
    xi = reshape(x{i},[R(i), N(i), R(i+1)]);
    xi = reshape(permute(xi, [2 1 3]),N(i),[]);
    temp = A{i}*xi;
    temp = reshape(temp,n_samples,R(i),R(i+1));
    temp2 = zeros(n_samples,R(i+1));
    for i2 = 1:R(i+1)
        temp2(:,i2) = sum(W.*temp(:,:,i2),2);
    end
    W = temp2;
end
b = W;
end