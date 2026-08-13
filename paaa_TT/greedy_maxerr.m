function [idx,max_val] = greedy_maxerr(x,init_idx)
%debug: ensure y{i} and yi are vectors, ensure the evaluated yi matches the
%corresponding indices, ensure the maximum indices are correct

[d,m,r] = TTsizes(x);
y = cell(d,1);
y{1} = x{1}(init_idx(1),:);
for i = 2:d-1
    xi = reshape(x{i},[r(i), m(i), r(i+1)]);
    y{i} = y{i-1}*reshape(xi(:,init_idx(i),:),[r(i), r(i+1)]);
end

idx = init_idx;
xi =  reshape(x{d},[r(d), m(d)]);
yi = y{d-1}*xi;
[max_val, ik] = max(abs(yi));
idx(d) = ik;
y{d} = xi(:,ik);

for i = d-1:-1:2
    xi = reshape(x{i},[r(i), m(i), r(i+1)]);
    yir = reshape(xi,[],r(i+1))*y{i+1}; %size r(i)*m(i)
    yir = reshape(yir,[r(i), m(i)]);
    yi = y{i-1}*yir;
    [max_val, ik] = max(abs(yi));
    idx(i) = ik;
    y{i} = xi(:, ik);
end

yi = x{1}*y{2};
[max_val, ik] = max(abs(yi));
idx(1) = ik;

end