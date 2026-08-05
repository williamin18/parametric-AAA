function [Y] = TTfull(X)
%not in original codes, written by William for testing

[N,I,R] = TTsizes(X);
n_vec = prod(I);
I_count = ones(N,1);
Y = zeros(I');
for j = 1:n_vec
    yj = 1;
    for m = 1:N
        yj = yj*X{m}(1+(I_count(m)-1)*R(m):I_count(m)*R(m),:);
    end
    Y(j) = yj;
    
    for k = 1:N
        if I_count(k) >= I(k)
            I_count(k) = 1;
        else
            I_count(k) = I_count(k)+1;
            break;
        end
    end
end