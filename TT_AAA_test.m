d = 4;
m = [25 25 25 25];

% Generate test data
x = linspace(-5,5,25);
y = linspace(-5,5,25);
z = linspace(-5,5,25);
w = linspace(-5,5,25);

samples = zeros(25,25,25,25);
for i = 1:25
    for j = 1:25
        for k = 1:25
            for p = 1:25
                samples(i,j,k,p) = (x(i) + y(j) + z(k) + w(p)) ./ (8 + cos(x(i)) + cos(y(j)) + cos(z(k)) + cos(w(p)));
            end
        end
    end
end



N = 2000;
idx = lhsdesign(N,4);
i1 = floor(idx(:,1)*25) + 1;
i2 = floor(idx(:,2)*25) + 1;
i3 = floor(idx(:,3)*25) + 1;
i4 = floor(idx(:,3)*25) + 1;

multi_idx = [i1,i2,i3,i4];
multi_idx = unique(multi_idx,'rows');
multi_idx = multi_idx(1:1000,:);

vec_idx =  sub2ind([25,25,25,25], ...
    multi_idx(:,1), ...
    multi_idx(:,2), ...
    multi_idx(:,3),...
    multi_idx(:,4));

sub_samples = samples(vec_idx);

r1_init = cell(d,1);
ref = [1 1 1 1];
for i = 1:d
    init_idx = repmat(ref,m(i),1);
    init_idx(:,i) = (1:m(i))';
    init_idx_cell = num2cell(init_idx,1);
    init_idx = sub2ind(m, init_idx_cell{:});
    
    r1_init{i} = samples(init_idx);
    [~,max_idx] = max(abs(r1_init{i} ));
    r1_init{i} = r1_init{i}/r1_init{i}(max_idx);
    ref(i) = r1_init{i};
end
