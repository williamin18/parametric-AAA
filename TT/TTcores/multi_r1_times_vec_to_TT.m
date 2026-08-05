function x = multi_r1_times_vec_to_TT(A,b)
%A'b = x
d = length(A);
[n_samples,m] = size(A{1});

x = cell(d,1);
x{1} = zeros(m,n_samples);
for j = 1:m
    x{1}(j,:) = (A{1}(:,j)').*(b.');
end
for i = 2:d-1
    [~,m] = size(A{i});
    x{i} = zeros(n_samples*m,n_samples);
    for j = 1:m
        x{i}((j-1)*n_samples+1:j*n_samples,:) = diag(A{i}(:,j)');
    end
end
[~,m] = size(A{d});
x{d} =  reshape(A{d},n_samples*m,1);
end