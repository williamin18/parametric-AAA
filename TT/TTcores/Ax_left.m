function [ylk,yl] = Ax_left(A,x,k)
%multiply rank 1 samples A with TT x for cores indices less than k, 
% yl{k} is a n*r_k matrix: n is number of samples, r_k is the rank of k's TT-core
[n_samples,~] = size(A{1});
[d,m,r] = TTsizes(x);


yl = cell(d,1);
yl{1} = ones(n_samples,1);
for i = 1:k-1
    xi = reshape(x{i},[r(i), m(i), r(i+1)]);
    xi = reshape(permute(xi, [2 1 3]),m(i),[]);
    Axi = A{i}*xi;
    Axi = reshape(Axi,n_samples,r(i),r(i+1));

    temp = zeros(n_samples,r(i+1));
    for j = 1:r(i+1)
        temp(:,j) = sum(yl{i}.*Axi(:,:,j),2);
    end
    yl{i+1} = temp;

end
ylk = yl{k};
end

