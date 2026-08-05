function [selector] = index2selector(idx,max_idx)
%INDEX2SELECTOR Summary of this function goes here
%   Detailed explanation goes here

[n,d]=size(idx);
if isscalar(max_idx)
    max_idx = max_idx*ones(d,1);
end

selector = cell(d,1);
rows = 1:n;
for i = 1:d
    selector{i} = zeros(n,max_idx(i));
    selector_idx = sub2ind([n,max_idx(i)],rows,idx(:,i)');
    selector{i}(selector_idx) = 1;
end

end

