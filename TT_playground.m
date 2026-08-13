x = linspace(-5,5,25);
y = linspace(-5,5,25);
z = linspace(-5,5,25);
samples = zeros(25,25,25);
for i = 1:25
    for j = 1:25
        for k = 1:25
            samples(i,j,k) = (x(i) + y(j) + z(k)) ./ (6 + cos(x(i)) + cos(y(j)) + cos(z(k)));
        end
    end
end

% Approximate data via low-rank p-AAA with error tolerance 1e-3 and tensor rank 3
lr_bf = paaa_TT(samples,{x,y,z},1e-3,3);
