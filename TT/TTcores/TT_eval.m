function [y] = TT_eval(x,I)
%TT_EVAL Summary of this function goes here
%   Detailed explanation goes here
[n_samples,d2] = size(I);
if d2 == 1
    I = I';
    n_samples = d2;
end


[d,N,R] = TTsizes(x);

y = x{1}(I(:,1),:);
for i = 2:d
    xi = reshape(x{i},[R(i), N(i), R(i+1)]);
    xi = reshape(permute(xi, [2 1 3]),N(i),[]);
    temp = xi(I(:,i),:);
    temp = reshape(temp,n_samples,R(i),R(i+1));
    temp2 = zeros(n_samples,R(i+1));
    for i2 = 1:R(i+1)
        temp2(:,i2) = sum(y.*temp(:,:,i2),2);
    end
    y = temp2;
end
end

